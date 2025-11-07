#!/usr/bin/env python
# coding: utf-8
"""
Spindle Positioning Simulation - Refactored Version
====================================================

Simulates spindle positioning during cell division across multiple cell types:
- Follicular epithelial cells (Drosophila)
- Neuroblasts (Drosophila)
- C. elegans embryos (PNC and spindle modes)
- Zebrafish endothelial cells

Author: Alikhan Yeltokov
Date: 2025
"""

from dataclasses import dataclass, field
from typing import Optional, Literal, Dict, Any, Tuple, List
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl
import cv2
import json
import random
import math
import os
from numpy import linalg as LA
from scipy.stats import norm
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from shapely.geometry import LineString, Polygon, Point
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors

# ============================================================================
# CONFIGURATION CLASSES
# ============================================================================

@dataclass
class SimulationConfig:
    """Core simulation parameters shared by all cell types."""
    # Force units are in pN, distances is in custom model units where 1 model unit = 10 µm (will change to micrometers in the next version), time in s
    # Time & discretization
    time_step: float
    total_time: float
    seed: int = 42
    
    # Spindle geometry
    spindle_half_length: float = 0.4
    spindle_width: float = 0.32
    spindle_angle_init: float = 0.0
    
    # Microtubules
    n_astral_mts: int = 100
    mt_mean_length: float = 1.0
    mt_length_stdev: float = 1.0
    
    # Forces & dynamics
    pull_force: float = 5.0
    push_force: float = 5.0
    growth_rate: float = 0.013
    shrink_rate: float = -0.027
    catastrophe_rate: float = 0.021
    rescue_rate: float = 0.029
    bind_rate: float = 0.03
    unbind_rate: float = 0.02
    
    # Physical constants
    viscosity: float = 100.0
    bending_rigidity: float = 0.2
    min_cortex_distance: float = 0.01
    push_distance: float = 0.0005
    max_interaction_distance: float = 0.02
    mt_min_length: float = 0.01
    repulsion_force: float = 0.0
    
    # Discretization
    n_cell_sides: int = 480
    mt_discretization: int = 2 #serves no function in the current version
    
    # Paths
    data_dir: Path = Path("./data")
    output_dir: Path = Path("./output")
    
    # Visualization
    show_force_vectors: bool = False
    save_plots: bool = True
    plot_interval: int = 20 # in time steps
    
    # Advanced
    length_distribution: Literal["gamma", "constant"] = "gamma"
    initial_state_distribution: Literal["random", "growing"] = "random"
    max_step: float = 0.1
    
    def get_probabilities(self) -> Tuple[float, float, float, float]:
        """Calculate transition probabilities from rates."""
        prob_catastrophe = 1 - np.exp(-self.catastrophe_rate * self.time_step)
        prob_rescue = 1 - np.exp(-self.rescue_rate * self.time_step)
        
        if self.pull_force > 0:
            prob_bind = 1 - np.exp(-self.bind_rate * self.time_step)
            prob_unbind = 1 - np.exp(-self.unbind_rate * self.time_step)
        else:
            prob_bind, prob_unbind = 0.0, 1.0
        
        return prob_catastrophe, prob_rescue, prob_bind, prob_unbind


@dataclass
class FollicleEpithelialConfig(SimulationConfig):
    cell_type: str = "follicle_epithelial"
    cell_radius_a: float = 0.5
    cell_radius_b: float = 0.5
    fg_density: int = 100
    astral_spread: float = 3 * np.pi / 2


@dataclass
class NeuroblastConfig(SimulationConfig):
    cell_type: str = "neuroblast"
    cell_radius_a: float = 0.5
    cell_radius_b: float = 0.5
    fg_density_basal: int = 100
    fg_density_apical: int = 100
    astral_spread: float = 3 * np.pi / 2


@dataclass
class CElegansPNCConfig(SimulationConfig):
    cell_type: str = "celegans_pnc"
    cell_radius_a: float = 2.5
    cell_radius_b: float = 1.5
    superellipse_n: float = 2.2
    fg_density_anterior: int = 40
    fg_density_posterior: int = 60
    astral_spread: float = np.pi
    catastrophe_rate: float = 0.014  # Override parent
    rescue_rate: float = 0.044        # Override parent
    spindle_angle_init: float = np.pi / 2  # Override parent


@dataclass
class CElegansSpindleConfig(SimulationConfig):
    cell_type: str = "celegans_spindle"
    cell_radius_a: float = 2.5
    cell_radius_b: float = 1.5
    superellipse_n: float = 2.2
    fg_density_anterior: int = 40
    fg_density_posterior: int = 60
    astral_spread: float = 3 * np.pi / 2
    catastrophe_rate: float = 0.014  # Override parent
    rescue_rate: float = 0.044        # Override parent


@dataclass
class ZebrafishEndoConfig(SimulationConfig):
    """Zebrafish endothelial cell from live imaging."""

    cell_type: str = "zebrafish_endo"
    endo_index: int = 1
    fg_distribution: Literal["uniform", "junctions"] = "junctions"
    fg_density: int = 10
    
    cell_radius_a: float = 1.0
    cell_radius_b: float = 1.0
    astral_spread: float = 1.5 * np.pi
    
    # Override parent's required field to be optional
    # If None, loads from Excel; if set, uses user value
    total_time: Optional[float] = None
    
    movie_info_file: Path = Path("./data/Movie_info.xlsx")
    cell_images_dir: Path = Path("./data/cells")
    spindle_images_dir: Path = Path("./data/spindles")
    junction_data_dir: Optional[Path] = None
    
    extras: Optional[Dict[str, Any]] = None
    
    def __post_init__(self):
        # Validation
        if self.fg_distribution == "junctions" and self.junction_data_dir is None:
            raise ValueError("junction_data_dir required for junctions mode")
        
        # Load parameters from Excel
        self._load_from_excel()
        
        # Only override total_time if user didn't provide it
        if self.total_time is None:
            self.total_time = ((self.extras["anaphase_on"] - self.extras["metaphase"]) 
                              * self.extras["frame_rate"])
    
    def _load_from_excel(self):
        """Load cell parameters from Movie_info.xlsx."""
        if not self.movie_info_file.exists():
            raise FileNotFoundError(f"Movie info not found: {self.movie_info_file}")
        
        df = pd.read_excel(self.movie_info_file, sheet_name="Alikhan", engine="openpyxl")
        if self.endo_index < 1 or self.endo_index > len(df):
            raise IndexError(f"endo_index {self.endo_index} out of range [1, {len(df)}]")
        
        row = df.iloc[self.endo_index - 1]
        
        self.extras = {
            "frame_rate": float(row["Frame rate"]),
            "anaphase_on": float(row["anaphase onset"]),
            "metaphase": float(row["metaphase"]),
            "rescale": float(row["model scale factor"]),
            "initial_angle": float(row["Initial angle"]),
            "final_angle": float(row["Final angle"]),
            "spindle_length": float(row["Spindle length"]),
            "spindle_width": float(row["Spindle width"]),
        }
        
        # Always override these from Excel (not optional)
        self.spindle_half_length = 0.1 * self.extras["spindle_length"] / 2.0
        self.spindle_width = 0.1 * self.extras["spindle_width"] / 2.0
        self.spindle_angle_init = np.deg2rad(self.extras["initial_angle"])


CellConfig = (FollicleEpithelialConfig | NeuroblastConfig | 
              CElegansPNCConfig | CElegansSpindleConfig | 
              ZebrafishEndoConfig)


# ============================================================================
# SIMULATION STATE
# ============================================================================

@dataclass
class SimulationState:
    """Complete state of the simulation at a given time."""
    
    time: float
    cell_boundary: np.ndarray
    fgs: np.ndarray  # Force generator positions
    spindle_poles: np.ndarray
    spindle_angle: float
    astral_mts: np.ndarray
    astral_angles: np.ndarray
    mt_states: np.ndarray  # +1 growing, -1 shrinking
    mt_lengths: np.ndarray
    which_bind: np.ndarray
    which_push: np.ndarray
    free_fgs: np.ndarray
    astral_which_fg: np.ndarray
    velocity_com: Optional[np.ndarray] = None


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def sgn(x):
    """Sign function."""
    return np.where(x > 0, 1, -1)


def gaussian_angles(
    n_points: int,
    mu: float = 0,
    sigma: float = 30,
    limits: Optional[Tuple[float, float]] = None
) -> np.ndarray:
    """Generate angles with Gaussian spacing."""
    uniform = np.linspace(0.001, 0.999, n_points)
    angles = norm.ppf(uniform, loc=mu, scale=sigma)
    
    if limits is not None:
        angles = angles[(angles >= limits[0]) & (angles <= limits[1])]
    
    return np.sort(angles)


def superellipse_points(
    angles_deg: np.ndarray,
    a: float,
    b: float,
    na: float
) -> np.ndarray:
    """Generate points on superellipse."""
    angles = np.deg2rad(angles_deg)
    x = (np.abs(np.cos(angles)) ** na) * a * sgn(np.cos(angles))
    y = (np.abs(np.sin(angles)) ** na) * b * sgn(np.sin(angles))
    return np.column_stack([x, y])


def point_in_polygon(point: np.ndarray, polygon: np.ndarray) -> bool:
    """Check if point is inside polygon."""
    x, y = point
    n = len(polygon)
    inside = False
    
    p1x, p1y = polygon[0]
    for i in range(n + 1):
        p2x, p2y = polygon[i % n]
        if y > min(p1y, p2y) and y <= max(p1y, p2y) and x <= max(p1x, p2x):
            if p1y != p2y:
                x_inters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                if p1x == p2x or x <= x_inters:
                    inside = not inside
        p1x, p1y = p2x, p2y
    
    return inside


def distance_matrix(points1: np.ndarray, points2: np.ndarray) -> np.ndarray:
    """Calculate pairwise distances."""
    d = points1[:, None, :] - points2[None, :, :]
    return np.linalg.norm(d, axis=2)


def calculate_perimeter(coordinates: np.ndarray) -> float:
    """Calculate perimeter of closed polygon."""
    distances = np.sqrt(np.sum(np.diff(coordinates, axis=0)**2, axis=1))
    closing = np.sqrt(np.sum((coordinates[0] - coordinates[-1])**2))
    return np.sum(distances) + closing


def calculate_partial_perimeter(
    coordinates: np.ndarray,
    i: int,
    j: int
) -> float:
    """Calculate perimeter between two indices."""
    if i > j:
        i, j = j, i
    subarray = coordinates[i:j+1]
    distances = np.sqrt(np.sum(np.diff(subarray, axis=0)**2, axis=1))
    closing = np.sqrt(np.sum((subarray[0] - subarray[-1])**2))
    return np.sum(distances) + closing


def slice_array(
    array: np.ndarray,
    start: int,
    end: int,
    num_points: int
) -> np.ndarray:
    """Extract evenly spaced points between indices."""
    if num_points < 2:
        raise ValueError("num_points must be >= 2")
    indices = np.linspace(start, end, num_points).astype(int)
    return array[indices]


def restructure_mt(mt: np.ndarray) -> np.ndarray:
    """Restructure MT to have evenly spaced points."""
    n_points = len(mt)
    restructured = np.ones_like(mt)
    restructured[:, 0] = np.linspace(mt[0, 0], mt[-1, 0], n_points)
    restructured[:, 1] = np.linspace(mt[0, 1], mt[-1, 1], n_points)
    return restructured


# ============================================================================
# IMAGE PROCESSING FUNCTIONS (for Zebrafish)
# ============================================================================

def get_contours(image_path: str) -> np.ndarray:
    """Extract contours from microscopy image using OpenCV."""
    image = cv2.imread(image_path)
    if image is None:
        raise FileNotFoundError(f"Could not load image: {image_path}")
    
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)
    edges = cv2.Canny(blurred, 50, 150)
    contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    smoothed_contours = []
    for contour in contours:
        epsilon = 0.0001 * cv2.arcLength(contour, True)
        smoothed_contour = cv2.approxPolyDP(contour, epsilon, True)
        smoothed_contours.append(smoothed_contour)
    
    all_points = np.concatenate([contour.squeeze() for contour in smoothed_contours])
    all_points[:, 1] = -all_points[:, 1]  # Flip y-axis
    
    return all_points


def enhance_cell(cell: np.ndarray, max_spacing: float = 0.005) -> np.ndarray:
    """Equalize spacing between cell contour vertices."""
    new_cell = []
    
    for i in range(len(cell)):
        next_i = 0 if i == len(cell) - 1 else i + 1
        dist = LA.norm(cell[next_i] - cell[i])
        
        if dist > max_spacing:
            new_cell.append(cell[i])
            n_points = int(dist / max_spacing)
            spread = np.linspace(0, dist, n_points)[1:-1]
            vec = (cell[next_i] - cell[i]) / dist
            
            for j in range(len(spread)):
                new_add = spread[j] * vec + cell[i]
                new_cell.append(new_add)
        else:
            new_cell.append(cell[i])
    
    return np.array(new_cell)


def smoothen_cell(cell: np.ndarray, sigma: float = 10) -> np.ndarray:
    """Apply Gaussian smoothing to cell contour."""
    x = cell[:, 0]
    y = cell[:, 1]
    
    x_smooth = gaussian_filter1d(x, sigma=sigma)
    y_smooth = gaussian_filter1d(y, sigma=sigma)
    
    return np.column_stack((x_smooth, y_smooth))


def close_contours(contours: np.ndarray, min_spacing: float = 0.005) -> np.ndarray:
    """Close gaps in contour by adding interpolated points."""
    start_point = contours[-1]
    end_point = contours[0]
    
    distance = np.linalg.norm(end_point - start_point)
    num_points_to_add = int(np.ceil(distance / min_spacing))
    
    if num_points_to_add > 0:
        x_new = np.linspace(start_point[0], end_point[0], num_points_to_add + 2)
        y_new = np.linspace(start_point[1], end_point[1], num_points_to_add + 2)
        
        x_new = x_new[1:-1]
        y_new = y_new[1:-1]
        
        new_points = np.column_stack((x_new, y_new))
        closed_contours = np.vstack([contours, new_points])
    else:
        closed_contours = contours
    
    return closed_contours


def resample_cell(traj: np.ndarray, n_points: int = 2400) -> np.ndarray:
    """Resample cell contour to have uniform spacing."""
    distances = np.cumsum(np.sqrt(np.sum(np.diff(traj, axis=0)**2, axis=1)))
    distances = np.insert(distances, 0, 0)
    
    fx = interp1d(distances, traj[:, 0], kind='linear')
    fy = interp1d(distances, traj[:, 1], kind='linear')
    
    new_distances = np.linspace(0, distances[-1], n_points)
    new_x = fx(new_distances)
    new_y = fy(new_distances)
    
    contour = np.column_stack((new_x, new_y))
    
    # Rotate so first point is at rightmost position
    target_point = np.array([np.max(contour[:, 0]), 0])
    distances_to_target = np.linalg.norm(contour - target_point, axis=1)
    closest_idx = np.argmin(distances_to_target)
    shifted_contour = np.roll(contour, -closest_idx, axis=0)
    
    return shifted_contour


def interpolate_contours(
    contour1: np.ndarray,
    contour2: np.ndarray,
    steps: int = 5
) -> List[np.ndarray]:
    """Generate interpolated contours between two shapes."""
    alphas = np.linspace(0, 1, steps)
    return [contour1 * (1 - alpha) + contour2 * alpha for alpha in alphas]


def transfer_free_fgs(new_length: int, old_array: np.ndarray) -> np.ndarray:
    """Transfer FG availability status to new array size."""
    new_array = np.zeros(new_length)
    min_length = min(new_length, len(old_array))
    new_array[:min_length] = old_array[:min_length]
    return new_array


def transfer_astral_which_fg(new_length: int, old_array: np.ndarray) -> np.ndarray:
    """Transfer MT-FG binding information to new array size."""
    new_array = np.zeros((new_length, 2))
    min_length = min(new_length, len(old_array))
    new_array[:min_length] = old_array[:min_length]
    return new_array


# ============================================================================
# GEOMETRY: CELL-CORTEX INTERSECTIONS
# ============================================================================

def intersect_cell(
    angle: float,
    start_point: np.ndarray,
    cell: np.ndarray
) -> Tuple[np.ndarray, str]:
    """
    Find intersection between MT ray and cell cortex using Shapely.
    
    Args:
        angle: MT angle (radians)
        start_point: MT origin (spindle pole)
        cell: Cell boundary coordinates (N, 2)
    
    Returns:
        intersection_point: (x, y) coordinates
        geom_type: 'Point', 'MultiPoint', or 'LineString'
    """
    # Create ray extending far beyond cell
    end_point = start_point + 4 * np.array([np.cos(angle), np.sin(angle)])
    
    mt_line = LineString([start_point, end_point])
    cell_line = LineString(cell)
    
    intersection = cell_line.intersection(mt_line)
    
    if intersection.geom_type == 'Point':
        return np.array([intersection.x, intersection.y]), intersection.geom_type
    
    elif intersection.geom_type == 'MultiPoint':
        # Multiple intersections - take closest to start
        points = list(intersection.geoms)
        distances = [LA.norm([start_point[0] - p.x, start_point[1] - p.y]) 
                    for p in points]
        closest = points[np.argmin(distances)]
        return np.array([closest.x, closest.y]), intersection.geom_type
    
    else:
        # Shouldn't happen, but fallback
        return end_point, intersection.geom_type


# ============================================================================
# FORCE GENERATOR (FG) POSITIONING
# ============================================================================

def make_fgs(
    config: CellConfig,
    cell: np.ndarray,
    frame: int = 1
) -> np.ndarray:
    """Generate force generator (motor protein) binding positions."""
    
    if isinstance(config, FollicleEpithelialConfig):
        return make_fgs_follicle(config)
    elif isinstance(config, NeuroblastConfig):
        return make_fgs_neuroblast(config)
    elif isinstance(config, CElegansPNCConfig):
        return make_fgs_celegans_pnc(config)
    elif isinstance(config, CElegansSpindleConfig):
        return make_fgs_celegans_spindle(config)
    elif isinstance(config, ZebrafishEndoConfig):
        return make_fgs_zebrafish(config, cell, frame)
    else:
        raise TypeError(f"Unknown config type: {type(config)}")


def make_fgs_follicle(config: FollicleEpithelialConfig) -> np.ndarray:
    """Basal/lateral FGs only."""
    n = config.fg_density
    a, b = config.cell_radius_a, config.cell_radius_b
    
    angles = gaussian_angles(n // 2, mu=10.7, sigma=30.3, limits=(-60, 60))
    fgs_right = np.column_stack([
        a * np.cos(np.deg2rad(angles)),
        b * np.sin(np.deg2rad(angles))
    ])
    
    fgs_left = fgs_right.copy()
    fgs_left[:, 0] *= -1
    
    return np.vstack([fgs_right[:-1], fgs_left])


def make_fgs_neuroblast(config: NeuroblastConfig) -> np.ndarray:
    """Apical FGs."""
    a, b = config.cell_radius_a, config.cell_radius_b
    
    # Apical
    n_apical = config.fg_density_apical
    angles_apical = gaussian_angles(n_apical, mu=90, sigma=30.3, limits=(45, 135))
    fgs_apical = np.column_stack([
        a * np.cos(np.deg2rad(angles_apical)),
        b * np.sin(np.deg2rad(angles_apical))
    ])
    
    return fgs_apical


def make_fgs_celegans_pnc(config: CElegansPNCConfig) -> np.ndarray:
    """Posterior-enriched distribution."""
    a, b, n = config.cell_radius_a, config.cell_radius_b, config.superellipse_n
    na = 2 / n
    
    n_ant = config.fg_density_anterior
    angles_ant = gaussian_angles(n_ant, mu=180, sigma=60, limits=(85, 275))
    fgs_ant = superellipse_points(angles_ant, a, b, na)
    
    n_post = config.fg_density_posterior
    angles_post = gaussian_angles(n_post, mu=0, sigma=60, limits=(-85, 85))
    fgs_post = superellipse_points(angles_post, a, b, na)
    
    return np.vstack([fgs_ant, fgs_post])


def make_fgs_celegans_spindle(config: CElegansSpindleConfig) -> np.ndarray:
    """Uniform distribution."""
    a, b, n = config.cell_radius_a, config.cell_radius_b, config.superellipse_n
    na = 2 / n
    
    n_ant = config.fg_density_anterior
    angles_ant = gaussian_angles(n_ant, mu=180, sigma=60, limits=(90, 270))
    fgs_ant = superellipse_points(angles_ant, a, b, na)
    
    n_post = config.fg_density_posterior
    angles_post = gaussian_angles(n_post, mu=0, sigma=60, limits=(-90, 90))
    fgs_post = superellipse_points(angles_post, a, b, na)
    
    return np.vstack([fgs_ant, fgs_post])


def make_fgs_zebrafish(
    config: ZebrafishEndoConfig,
    cell: np.ndarray,
    frame: int
) -> np.ndarray:
    """Uniform or junction-localized FGs."""
    
    if config.fg_distribution == "uniform":
        perimeter = calculate_perimeter(cell)
        n_fgs = int(config.fg_density * perimeter)
        spacing = max(1, len(cell) // n_fgs)
        return cell[::spacing]
    
    else:  # junctions
        juncs = load_junctions(config)
        junc_frame = juncs[frame - 1]
        
        dist1 = distance_matrix(np.array([junc_frame[1]]), cell)
        end_idx = np.argmin(dist1[0])
        
        dist2 = distance_matrix(np.array([junc_frame[0]]), cell)
        start_idx = np.argmin(dist2[0])
        
        partial_perim = calculate_partial_perimeter(cell, start_idx, end_idx)
        n_fgs = int(config.fg_density * partial_perim)
        
        return slice_array(cell, start_idx, end_idx, n_fgs)


def load_junctions(config: ZebrafishEndoConfig) -> np.ndarray:
    """Load junction coordinates for zebrafish cells."""
    junction_file = (config.junction_data_dir / 
                    f"C{config.endo_index}_Mask_Movie_Junctions.txt")
    
    if not junction_file.exists():
        raise FileNotFoundError(f"Junction data not found: {junction_file}")
    
    data = np.loadtxt(junction_file)
    data[:, 1] = -data[:, 1]
    return data.reshape(-1, 2, 2)


# ============================================================================
# CELL GEOMETRY
# ============================================================================

def get_real_cell_zebrafish(
    config: ZebrafishEndoConfig,
    t_time: float = 0
) -> np.ndarray:
    """Load and process cell contours from microscopy images."""
    frame = int(config.extras["metaphase"] + t_time / config.extras["frame_rate"])
    
    # Load cell image
    image_path = config.cell_images_dir / f"cell_{config.endo_index}" / f"Mask_{frame}.jpg"
    if not image_path.exists():
        raise FileNotFoundError(f"Cell image not found: {image_path}")
    
    # Extract and process contours
    all_points = get_contours(str(image_path))
    rescale = config.extras["rescale"]
    
    cell_input = np.zeros((len(all_points), 2))
    cell_input[:, 0] = all_points[:, 0] / rescale
    cell_input[:, 1] = all_points[:, 1] / rescale
    
    cell = enhance_cell(cell_input)
    
    # Get spindle center for recentering (0,0 at spindle center at time 0)
    spindle_image_dir = config.spindle_images_dir / f"spindle_{config.endo_index}"
    if config.endo_index > 9:
        spindle_path = spindle_image_dir / f"Mask_{int(config.extras['metaphase'])}.jpg"
    else:
        spindle_path = spindle_image_dir / f"Spindle_{int(config.extras['metaphase'])}.jpg"
    
    if not spindle_path.exists():
        raise FileNotFoundError(f"Spindle image not found: {spindle_path}")
    
    spindle_contour = get_contours(str(spindle_path))
    spindle_input = np.zeros((len(spindle_contour), 2))
    spindle_input[:, 0] = spindle_contour[:, 0] / rescale
    spindle_input[:, 1] = spindle_contour[:, 1] / rescale
    
    shapely_string = LineString(spindle_input)
    spindle_center = np.array([shapely_string.centroid.x, shapely_string.centroid.y])
    
    # Center cell on spindle
    cell = cell - spindle_center
    cell = smoothen_cell(cell)
    cell = close_contours(cell)
    cell = resample_cell(cell, n_points=2400)
    
    return cell


def update_cell_zebrafish(config: ZebrafishEndoConfig, t_time: float) -> np.ndarray:
    """Update cell shape with interpolation between frames."""
    frame_rate = config.extras["frame_rate"]
    steps = int(frame_rate)
    
    cell_before = get_real_cell_zebrafish(config, t_time)
    cell_after = get_real_cell_zebrafish(config, t_time + frame_rate)
    
    interpolated = interpolate_contours(cell_before, cell_after, steps)
    current = int((t_time / frame_rate % 1) * steps)
    
    return interpolated[current]


def make_cell_boundary(config: CellConfig, t_time: float = 0) -> np.ndarray:
    """Generate cell boundary coordinates."""
    
    theta = np.linspace(0, 2*np.pi, config.n_cell_sides + 1)[:-1]
    
    if isinstance(config, ZebrafishEndoConfig):
        return get_real_cell_zebrafish(config, t_time)
    
    elif isinstance(config, (CElegansPNCConfig, CElegansSpindleConfig)):
        a, b = config.cell_radius_a, config.cell_radius_b
        n = config.superellipse_n
        na = 2 / n
        cell = np.zeros((config.n_cell_sides, 2))
        cell[:, 0] = (np.abs(np.cos(theta)) ** na) * a * sgn(np.cos(theta))
        cell[:, 1] = (np.abs(np.sin(theta)) ** na) * b * sgn(np.sin(theta))
        return cell
    
    else:  # Follicle or Neuroblast
        cell = np.zeros((config.n_cell_sides, 2))
        cell[:, 0] = config.cell_radius_a * np.cos(theta)
        cell[:, 1] = config.cell_radius_b * np.sin(theta)
        return cell


def initialize_spindle(
    config: CellConfig,
    cell: np.ndarray
) -> Tuple[np.ndarray, float]:
    """
    Initialize spindle pole positions and angle.
    
    Returns:
        spindle_poles: (2, 2) array of pole positions
        spindle_angle: Adjusted angle (radians)
    """
    spindle_poles = np.zeros((2, 2))
    spindle_angle = config.spindle_angle_init
    r = config.spindle_half_length
    
    # Cell-type specific displacement
    displacement = np.array([0.0, 0.0])
    if isinstance(config, CElegansPNCConfig):
        displacement = np.array([1.0, 0.0]) # Shift toward posterior
    
    spindle_poles[0] = [
        r * np.cos(spindle_angle) + displacement[0],
        r * np.sin(spindle_angle) + displacement[1]
    ]
    spindle_poles[1] = [
        r * np.cos(spindle_angle + np.pi) + displacement[0],
        r * np.sin(spindle_angle + np.pi) + displacement[1]
    ]
    
    # Adjust angle if poles outside cell
    attempts = 0
    while (not point_in_polygon(spindle_poles[0], cell) or 
           not point_in_polygon(spindle_poles[1], cell)):
        spindle_angle += np.pi / 12
        spindle_poles[0] = [
            r * np.cos(spindle_angle) + displacement[0],
            r * np.sin(spindle_angle) + displacement[1]
        ]
        spindle_poles[1] = [
            r * np.cos(spindle_angle + np.pi) + displacement[0],
            r * np.sin(spindle_angle + np.pi) + displacement[1]
        ]
        attempts += 1
        if attempts > 12:
            print("Warning: Could not fit spindle inside cell")
            break
    
    return spindle_poles, spindle_angle


# ============================================================================
# MICROTUBULE INITIALIZATION & DYNAMICS
# ============================================================================

def initialize_microtubules(
    config: CellConfig,
    cell: np.ndarray,
    spindle_poles: np.ndarray,
    spindle_angle: float,
    fgs: np.ndarray
) -> Tuple:
    """
    Initialize astral microtubules.
    
    Returns:
        astral_mts: (2, n_astro, discr, 2) MT coordinates
        astral_angles: (2, n_astro) MT angles
        mt_states: (2, n_astro) growth states (+1/-1)
        which_push: (2, n_astro) pushing status
        which_bind: (2, n_astro) binding status
        free_fgs: (n_fgs,) available FG flags
        astral_which_fg: (n_fgs, 2) FG assignments
        mt_lengths: (2, n_astro) current MT lengths
    """
    n_astro = config.n_astral_mts
    n_fgs = len(fgs)
    
    # Initialize arrays
    free_fgs = np.zeros(n_fgs)
    astral_which_fg = np.zeros((n_fgs, 2)) - 1
    which_bind = np.zeros((2, n_astro))
    which_push = np.zeros((2, n_astro))
    astral_mts = np.zeros((2, n_astro, config.mt_discretization, 2))
    
    # MT angles distributed around each pole
    astral_angles = np.zeros((2, n_astro))
    spread = config.astral_spread
    astral_angles[0] = np.linspace(
        spindle_angle - spread/2, 
        spindle_angle + spread/2, 
        n_astro
    )
    astral_angles[1] = np.linspace(
        spindle_angle + np.pi - spread/2,
        spindle_angle + np.pi + spread/2, 
        n_astro
    )
    
    # Initialize lengths
    if config.length_distribution == "gamma":
        mt_lengths = np.random.gamma(
            config.mt_mean_length, 
            config.mt_length_stdev, 
            size=(2, n_astro)
        )
    else:
        mt_lengths = config.mt_mean_length * np.ones((2, n_astro))
    
    # Initialize growth states
    prob_cat, prob_res, _, _ = config.get_probabilities()
    if config.initial_state_distribution == "random":
        p_grow = prob_res / (prob_res + prob_cat)
        p_shrink = prob_cat / (prob_res + prob_cat)
        mt_states = np.random.choice(
            [-1, 1], 
            size=(2, n_astro), 
            p=[p_shrink, p_grow]
        )
    else:
        mt_states = np.ones((2, n_astro))
    
    # Build individual MTs
    for i in range(2):
        for j in range(n_astro):
            end = grow_mt_end(
                config, 
                astral_angles[i, j], 
                spindle_poles[i], 
                cell, 
                mt_lengths[i, j]
            )
            
            astral_mts[i, j, :, 0] = np.linspace(
                spindle_poles[i, 0], 
                end[0], 
                config.mt_discretization
            )
            astral_mts[i, j, :, 1] = np.linspace(
                spindle_poles[i, 1], 
                end[1], 
                config.mt_discretization
            )
            
            # Check initial binding/pushing
            (which_bind[i, j], which_push[i, j], mt_states[i, j], 
             free_fgs, astral_which_fg) = check_mt_interactions(
                config, i, j, astral_mts[i, j, -1], spindle_poles,
                fgs, free_fgs, astral_which_fg, mt_states[i, j],
                astral_angles[i, j], cell
            )
            
            mt_lengths[i, j] = LA.norm(astral_mts[i, j, -1] - spindle_poles[i])
    
    return (astral_mts, astral_angles, mt_states, which_push, which_bind,
            free_fgs, astral_which_fg, mt_lengths)


def grow_mt_end(
    config: CellConfig,
    angle: float,
    pole: np.ndarray,
    cell: np.ndarray,
    length: float
) -> np.ndarray:
    """
    Calculate MT end position after growth.
    
    Args:
        config: Simulation configuration
        angle: MT angle (radians)
        pole: Spindle pole position
        cell: Cell boundary
        length: Target MT length
    
    Returns:
        end_point: MT tip position (x, y)
    """
    # Find cortex intersection
    intersect, _ = intersect_cell(angle, pole, cell)
    
    # Proposed end after growth
    end = pole + (length + config.growth_rate * config.time_step) * np.array([
        np.cos(angle), np.sin(angle)
    ])
    
    # Check if end exceeds cortex
    len_to_intersect = LA.norm(intersect - pole)
    len_to_end = LA.norm(end - pole)
    
    if len_to_end < len_to_intersect:
        return end
    else:
        return intersect


def check_mt_interactions(
    config: CellConfig,
    pole_idx: int,
    mt_idx: int,
    mt_tip: np.ndarray,
    spindle_poles: np.ndarray,
    fgs: np.ndarray,
    free_fgs: np.ndarray,
    astral_which_fg: np.ndarray,
    state: int,
    angle: float,
    cell: np.ndarray
) -> Tuple[int, int, int, np.ndarray, np.ndarray]:
    """
    Check if MT binds to FG or pushes against cortex.
    
    Returns:
        bind: 1 if bound, 0 otherwise
        push: 1 if pushing, 0 otherwise
        state: Updated state
        free_fgs: Updated availability
        astral_which_fg: Updated assignments
    """
    bind = 0
    push = 0
    
    prob_cat, prob_res, prob_bind, prob_unbind = config.get_probabilities()
    
    # Check binding to FGs
    dist = distance_matrix(np.array([mt_tip]), fgs)
    min_dist_idx = np.argmin(dist[0])
    min_dist = dist[0, min_dist_idx]
    
    if (min_dist <= config.max_interaction_distance and 
        free_fgs[min_dist_idx] == 0 and 
        random.uniform(0, 1) <= prob_bind):
        
        free_fgs[min_dist_idx] = 1
        bind = 1
        state = -1  # Switch to shrinking
        astral_which_fg[min_dist_idx, 0] = pole_idx
        astral_which_fg[min_dist_idx, 1] = mt_idx
    
    else:
        # Check pushing against cortex
        cortex_point, _ = intersect_cell(angle, spindle_poles[pole_idx], cell)
        dist_to_cortex = LA.norm(cortex_point - mt_tip)
        
        if dist_to_cortex <= config.push_distance:
            push = 1
            state = 1  # Keep growing
    
    return bind, push, state, free_fgs, astral_which_fg


# ============================================================================
# FORCE CALCULATIONS
# ============================================================================

def calculate_forces(
    config: CellConfig,
    state: SimulationState
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate push and pull forces on each MT.
    
    Returns:
        push_forces: (2, n_astro, 2) push force vectors
        pull_forces: (2, n_astro, 2) pull force vectors
    """
    n_astro = config.n_astral_mts
    push_forces = np.zeros((2, n_astro, 2))
    pull_forces = np.zeros((2, n_astro, 2))
    
    for i in range(2):
        for j in range(n_astro):
            # Unit vector along MT
            mt_vec = (state.astral_mts[i, j, -1] - state.spindle_poles[i])
            mt_length = LA.norm(mt_vec)
            
            if mt_length > 0:
                unit_vec = mt_vec / mt_length
            else:
                continue
            
            # Pulling force
            if state.which_bind[i, j] == 1 and state.which_push[i, j] == 0:
                pull_forces[i, j] = config.pull_force * unit_vec
            
            # Pushing force (with Euler buckling limit)
            elif state.which_bind[i, j] == 0 and state.which_push[i, j] == 1:
                buckling_force = (config.bending_rigidity * np.pi**2 / 
                                 (mt_length**2))
                actual_push = min(config.push_force, buckling_force)
                push_forces[i, j] = -actual_push * unit_vec
    
    return push_forces, pull_forces


def calculate_torque(
    config: CellConfig,
    state: SimulationState,
    force1: np.ndarray,
    force2: np.ndarray
) -> float:
    """
    Calculate net torque on spindle.
    
    Args:
        config: Configuration
        state: Current state
        force1: Net force on pole 1
        force2: Net force on pole 2
    
    Returns:
        torque: Net torque (scalar)
    """
    r = config.spindle_half_length
    
    # Spindle axis unit vector
    poles_vec = state.spindle_poles[0] - state.spindle_poles[1]
    if LA.norm(poles_vec) > 0:
        poles_unit = poles_vec / (2 * r)
    else:
        return 0.0
    
    # Transverse components
    force1_t = force1 - np.dot(force1, poles_unit) * poles_unit
    force2_t = force2 - np.dot(force2, -poles_unit) * (-poles_unit)
    
    # Signs from cross product
    sign1 = np.sign(np.cross(poles_unit, force1_t))
    sign2 = np.sign(np.cross(-poles_unit, force2_t))
    
    # Torque = r × F
    torque = sign1 * LA.norm(force1_t) * r + sign2 * LA.norm(force2_t) * r
    
    return torque


# ============================================================================
# SPINDLE MOTION
# ============================================================================

def update_spindle_position(
    config: CellConfig,
    state: SimulationState,
    push_forces: np.ndarray,
    pull_forces: np.ndarray
) -> Tuple[np.ndarray, float, np.ndarray]:
    """
    Update spindle position and angle based on forces.
    
    Returns:
        new_poles: Updated pole positions
        new_angle: Updated angle
        velocity: Center of mass velocity
    """
    # Net forces on each pole
    force1 = np.sum(pull_forces[0] + push_forces[0], axis=0)
    force2 = np.sum(pull_forces[1] + push_forces[1], axis=0)
    
    # Add repulsion from cortex if needed
    force1 += calculate_cortex_repulsion(config, state, pole_idx=0)
    force2 += calculate_cortex_repulsion(config, state, pole_idx=1)
    
    net_force = force1 + force2
    torque = calculate_torque(config, state, force1, force2)
    
    # Translational and rotational velocities (Stokes drag)
    r = config.spindle_half_length
    mu = config.viscosity
    
    V = net_force / (6 * np.pi * mu * r)  # Translation
    Omega = torque / (8 * np.pi * mu * r**3)  # Rotation
    
    # Proposed new configuration
    new_poles = state.spindle_poles + V * config.time_step
    new_angle = state.spindle_angle + Omega * config.time_step
    
    # Recompute poles based on new angle (maintain spindle length)
    com = np.mean(new_poles, axis=0)
    new_poles[0] = com + r * np.array([np.cos(new_angle), np.sin(new_angle)])
    new_poles[1] = com + r * np.array([np.cos(new_angle + np.pi), np.sin(new_angle + np.pi)])
    
    # Check if valid (inside cell and away from cortex)
    if not check_spindle_valid(config, new_poles, new_angle, state.cell_boundary):
        # If invalid, try sliding along boundary
        new_poles, new_angle, success = slide_spindle(
            config, state, new_poles, new_angle, V, Omega
        )
        if not success:
            # Keep old position
            new_poles = state.spindle_poles.copy()
            new_angle = state.spindle_angle
            V = np.zeros(2)
    
    return new_poles, new_angle, V


def calculate_cortex_repulsion(
    config: CellConfig,
    state: SimulationState,
    pole_idx: int
) -> np.ndarray:
    """Calculate repulsive force from nearby cortex."""
    pole = state.spindle_poles[pole_idx]
    
    dist = distance_matrix(np.array([pole]), state.cell_boundary)
    min_dist = np.min(dist)
    
    if min_dist <= 1.5 * config.min_cortex_distance:
        # Repulsion away from cell center
        repel_vec = -pole / LA.norm(pole) if LA.norm(pole) > 0 else np.zeros(2)
        return (config.repulsion_force / min_dist) * repel_vec
    
    return np.zeros(2)


def check_spindle_valid(
    config: CellConfig,
    poles: np.ndarray,
    angle: float,
    cell: np.ndarray
) -> bool:
    """Check if spindle configuration is valid."""
    # Both poles inside cell
    if not point_in_polygon(poles[0], cell):
        return False
    if not point_in_polygon(poles[1], cell):
        return False
    
    # Minimum distance from cortex
    spindle_body = generate_spindle_body(config, poles, angle)
    cell_poly = Polygon(cell)
    spindle_poly = Polygon(spindle_body)
    
    return cell_poly.buffer(-config.min_cortex_distance).contains(spindle_poly)


def generate_spindle_body(
    config: CellConfig,
    poles: np.ndarray,
    angle: float
) -> np.ndarray:
    """Generate spindle body outline as superellipse."""
    theta = np.linspace(0, 2*np.pi, 36, endpoint=False)
    com = np.mean(poles, axis=0)
    
    # Superellipse
    a, b = config.spindle_half_length, config.spindle_width
    n = 1.4
    na = 2 / n
    
    x = (np.abs(np.cos(theta)) ** na) * a * sgn(np.cos(theta))
    y = (np.abs(np.sin(theta)) ** na) * b * sgn(np.sin(theta))
    spindle = np.vstack((x, y)).T + com
    
    # Rotate
    R = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle),  np.cos(angle)]
    ])
    spindle = (spindle - com) @ R.T + com
    
    return spindle


def slide_spindle(
    config: CellConfig,
    state: SimulationState,
    target_poles: np.ndarray,
    target_angle: float,
    V: np.ndarray,
    Omega: float
) -> Tuple[np.ndarray, float, bool]:
    """
    Try to slide spindle along boundary if blocked.
    
    Returns:
        poles: Best achievable position
        angle: Best achievable angle
        success: True if any movement was possible
    """
    # Binary search for maximum translation
    low, high = 0, 1
    for _ in range(20):
        frac = (low + high) / 2
        test_poles = state.spindle_poles + frac * V * config.time_step
        
        if check_spindle_valid(config, test_poles, state.spindle_angle, state.cell_boundary):
            low = frac
        else:
            high = frac
    
    trans_poles = state.spindle_poles + low * V * config.time_step
    
    # Binary search for maximum rotation
    low_rot, high_rot = 0, 1
    for _ in range(20):
        frac = (low_rot + high_rot) / 2
        test_angle = state.spindle_angle + frac * Omega * config.time_step
        
        com = np.mean(trans_poles, axis=0)
        r = config.spindle_half_length
        test_poles_rot = np.array([
            com + r * np.array([np.cos(test_angle), np.sin(test_angle)]),
            com + r * np.array([np.cos(test_angle + np.pi), np.sin(test_angle + np.pi)])
        ])
        
        if check_spindle_valid(config, test_poles_rot, test_angle, state.cell_boundary):
            low_rot = frac
        else:
            high_rot = frac
    
    final_angle = state.spindle_angle + low_rot * Omega * config.time_step
    com = np.mean(trans_poles, axis=0)
    r = config.spindle_half_length
    final_poles = np.array([
        com + r * np.array([np.cos(final_angle), np.sin(final_angle)]),
        com + r * np.array([np.cos(final_angle + np.pi), np.sin(final_angle + np.pi)])
    ])
    
    moved = (low > 0.01 or low_rot > 0.01)
    
    return final_poles, final_angle, moved


# ============================================================================
# MICROTUBULE DYNAMICS UPDATE
# ============================================================================

def update_microtubules(
    config: CellConfig,
    state: SimulationState,
    delta_angle: float
) -> SimulationState:
    """
    Update all microtubule states, lengths, and interactions.
    
    Args:
        config: Configuration
        state: Current state
        delta_angle: Change in spindle angle this timestep
    
    Returns:
        Updated state
    """
    prob_cat, prob_res, prob_bind, prob_unbind = config.get_probabilities()
    
    # Update MT angles with spindle rotation
    new_astral_angles = state.astral_angles + delta_angle
    
    new_astral_mts = state.astral_mts.copy()
    new_mt_states = state.mt_states.copy()
    new_which_bind = state.which_bind.copy()
    new_which_push = state.which_push.copy()
    new_free_fgs = state.free_fgs.copy()
    new_astral_which_fg = state.astral_which_fg.copy()
    new_mt_lengths = state.mt_lengths.copy()
    
    for i in range(2):
        for j in range(config.n_astral_mts):
            # Update minus end to new pole position
            new_astral_mts[i, j, 0] = state.spindle_poles[i]
            
            # Handle bound MTs
            if new_which_bind[i, j] == 1:
                # Check for unbinding
                cortex_point, geom_type = intersect_cell(
                    new_astral_angles[i, j],
                    state.spindle_poles[i],
                    state.cell_boundary
                )
                
                if (random.uniform(0, 1) <= prob_unbind or 
                    geom_type == 'MultiPoint'):
                    # Unbind
                    new_which_bind[i, j] = 0
                    fg_idx = np.where(
                        (new_astral_which_fg[:, 0] == i) & 
                        (new_astral_which_fg[:, 1] == j)
                    )[0]
                    if len(fg_idx) > 0:
                        new_free_fgs[fg_idx[0]] = 0
                    new_mt_states[i, j] = -1
                    
                    if geom_type == 'MultiPoint':
                        new_astral_mts[i, j, -1] = cortex_point
                        new_astral_mts[i, j] = restructure_mt(new_astral_mts[i, j])
            
            # Handle free MTs
            else:
                # Catastrophe/rescue transitions
                rand_val = random.uniform(0, 1)
                
                if new_mt_states[i, j] > 0 and rand_val <= prob_cat:
                    new_mt_states[i, j] = -1
                    rate = config.shrink_rate
                elif new_mt_states[i, j] < 0 and rand_val <= prob_res:
                    new_mt_states[i, j] = 1
                    rate = config.growth_rate
                else:
                    rate = (config.shrink_rate if new_mt_states[i, j] == -1 
                           else config.growth_rate)
                
                # Update length
                if new_mt_states[i, j] == -1:
                    # Shrinking
                    old_tip = (state.spindle_poles[i] + 
                              new_mt_lengths[i, j] * np.array([
                                  np.cos(new_astral_angles[i, j]),
                                  np.sin(new_astral_angles[i, j])
                              ]))
                    new_tip = old_tip + rate * config.time_step * np.array([
                        np.cos(new_astral_angles[i, j]),
                        np.sin(new_astral_angles[i, j])
                    ])
                    
                    # Check doesn't exceed cortex
                    cortex_point, _ = intersect_cell(
                        new_astral_angles[i, j],
                        state.spindle_poles[i],
                        state.cell_boundary
                    )
                    
                    if LA.norm(new_tip - state.spindle_poles[i]) > LA.norm(cortex_point - state.spindle_poles[i]):
                        new_tip = cortex_point
                    
                    new_astral_mts[i, j, -1] = new_tip
                    new_astral_mts[i, j] = restructure_mt(new_astral_mts[i, j])
                    
                    # Check if too short
                    if LA.norm(new_tip - state.spindle_poles[i]) < config.mt_min_length:
                        # Nucleate new MT
                        new_length = (config.mt_mean_length if config.length_distribution == "constant"
                                     else np.random.gamma(config.mt_mean_length, config.mt_length_stdev))
                        new_tip = grow_mt_end(
                            config,
                            new_astral_angles[i, j],
                            state.spindle_poles[i],
                            state.cell_boundary,
                            new_length
                        )
                        new_astral_mts[i, j, -1] = new_tip
                        new_astral_mts[i, j] = restructure_mt(new_astral_mts[i, j])
                        new_mt_states[i, j] = 1
                
                else:
                    # Growing
                    new_tip = grow_mt_end(
                        config,
                        new_astral_angles[i, j],
                        state.spindle_poles[i],
                        state.cell_boundary,
                        new_mt_lengths[i, j]
                    )
                    new_astral_mts[i, j, -1] = new_tip
                    new_astral_mts[i, j] = restructure_mt(new_astral_mts[i, j])
                
                # Update length
                new_mt_lengths[i, j] = LA.norm(
                    new_astral_mts[i, j, -1] - state.spindle_poles[i]
                )
                
                # Check for new binding/pushing
                (new_which_bind[i, j], new_which_push[i, j], new_mt_states[i, j],
                 new_free_fgs, new_astral_which_fg) = check_mt_interactions(
                    config, i, j, new_astral_mts[i, j, -1], state.spindle_poles,
                    state.fgs, new_free_fgs, new_astral_which_fg, new_mt_states[i, j],
                    new_astral_angles[i, j], state.cell_boundary
                )
    
    # Create updated state
    new_state = SimulationState(
        time=state.time + config.time_step,
        cell_boundary=state.cell_boundary,
        fgs=state.fgs,
        spindle_poles=state.spindle_poles,
        spindle_angle=state.spindle_angle,
        astral_mts=new_astral_mts,
        astral_angles=new_astral_angles,
        mt_states=new_mt_states,
        mt_lengths=new_mt_lengths,
        which_bind=new_which_bind,
        which_push=new_which_push,
        free_fgs=new_free_fgs,
        astral_which_fg=new_astral_which_fg,
        velocity_com=state.velocity_com
    )
    
    return new_state


# ============================================================================
# VISUALIZATION
# ============================================================================

def plot_cell(
    config: CellConfig,
    state: SimulationState,
    push_forces: np.ndarray,
    pull_forces: np.ndarray,
    step: int,
    output_dir: Path
):
    """
    Plot cell, spindle, MTs, and forces (original visualization style).
    
    Args:
        config: Simulation configuration
        state: Current simulation state
        push_forces: Push force vectors (2, n_astro, 2)
        pull_forces: Pull force vectors (2, n_astro, 2)
        step: Current timestep number
        output_dir: Directory to save plots
    """
    if not config.save_plots:
        return
    
    # Calculate force statistics
    pull_total = np.sum(pull_forces[0] + pull_forces[1], axis=0)
    push_total = np.sum(push_forces[0] + push_forces[1], axis=0)
    
    count_pulling = np.sum(state.which_bind == 1, axis=1)
    count_pushing = np.sum(state.which_push == 1, axis=1)
    
    total_pull_mag = LA.norm(pull_total)
    total_push_mag = LA.norm(push_total)
    
    if (total_pull_mag + total_push_mag) == 0:
        ratio = 0
    else:
        ratio = 100 * total_pull_mag / (total_pull_mag + total_push_mag)
    
    # Set up plot boundaries based on cell type
    cell_type = config.cell_type
    a, b = config.cell_radius_a, config.cell_radius_b
    
    if cell_type in ['follicle_epithelial', 'neuroblast']:
        left = min(state.cell_boundary[:, 0]) - 1.5*a
        right = max(state.cell_boundary[:, 0]) + 1.5*a
        top = max(state.cell_boundary[:, 1]) + 0.5
        bottom = min(state.cell_boundary[:, 1]) - 0.5
    elif cell_type in ['celegans_pnc', 'celegans_spindle']:
        left = min(state.cell_boundary[:, 0]) - 0.5*a
        right = max(state.cell_boundary[:, 0]) + 0.5*a
        top = max(state.cell_boundary[:, 1]) + 1.0
        bottom = min(state.cell_boundary[:, 1]) - 1.0
    else:  # zebrafish_endo
        left = min(state.cell_boundary[:, 0]) - 2
        right = max(state.cell_boundary[:, 0]) + 2
        top = max(state.cell_boundary[:, 1]) + 1
        bottom = min(state.cell_boundary[:, 1]) - 1
    
    # Create figure
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.set_xlim([left, right])
    ax.set_ylim([bottom, top])
    
    # Don't show axes
    ax.axis('off')
    
    # Plot cell boundary
    ax.plot(state.cell_boundary[:, 0], state.cell_boundary[:, 1], 
            color='darkgrey', label=f'Time = {state.time:.1f} s', linewidth=5, zorder=10)
    ax.plot([state.cell_boundary[0, 0], state.cell_boundary[-1, 0]], 
            [state.cell_boundary[0, 1], state.cell_boundary[-1, 1]], 
            color='darkgrey', label=f'Angle = {np.rad2deg(state.spindle_angle):.2f}°', linewidth=5)
    
    # Plot spindle poles
    ax.scatter([state.spindle_poles[0, 0], state.spindle_poles[1, 0]],
               [state.spindle_poles[0, 1], state.spindle_poles[1, 1]],
               color='yellow', s=150, edgecolors='k', zorder=6)
    
    # Plot FGs (force generators)
    ax.scatter(state.fgs[:, 0], state.fgs[:, 1], 
               label=f'N_MTs = {config.n_astral_mts*2}, N_FGs = {len(state.fgs)}',
               color='salmon', s=55, edgecolors='dimgrey', marker="8", zorder=50)
    
    # Plot MTs with color coding
    for i in range(2):
        for j in range(config.n_astral_mts):
            mt = state.astral_mts[i, j]
            
            if state.which_bind[i, j] == 1:  # Bound (pulling)
                color = 'red'
                zorder = 5
            elif state.mt_states[i, j] == -1:  # Shrinking
                color = 'tab:cyan'
                zorder = 3
            elif state.which_push[i, j] == 1:  # Pushing
                color = 'darkgreen' if i == 0 else 'green'
                zorder = 4
            else:  # Growing (free)
                color = 'slateblue'
                zorder = 2
            
            ax.plot(mt[:, 0], mt[:, 1], color=color, zorder=zorder)
    
    # Calculate net forces on each pole (needed for both plotting and text)
    force1 = np.sum(pull_forces[0] + push_forces[0], axis=0)
    force2 = np.sum(pull_forces[1] + push_forces[1], axis=0)
    
    pull_1 = np.sum(pull_forces[0], axis=0)
    push_1 = np.sum(push_forces[0], axis=0)
    pull_2 = np.sum(pull_forces[1], axis=0)
    push_2 = np.sum(push_forces[1], axis=0)
    
    # Force vectors (if enabled)
    if config.show_force_vectors:
        kap = 0.01  # Scaling factor for force arrows
        
        # Individual forces
        ax.plot([state.spindle_poles[0, 0], state.spindle_poles[0, 0] - kap*push_1[0]],
                [state.spindle_poles[0, 1], state.spindle_poles[0, 1] - kap*push_1[1]],
                'tab:red', linewidth=4, label=f'F_push1 = {LA.norm(push_1):.2f}')
        ax.plot([state.spindle_poles[0, 0], state.spindle_poles[0, 0] + kap*pull_1[0]],
                [state.spindle_poles[0, 1], state.spindle_poles[0, 1] + kap*pull_1[1]],
                'tab:pink', linewidth=4, label=f'F_pull1 = {LA.norm(pull_1):.2f}')
        
        ax.plot([state.spindle_poles[1, 0], state.spindle_poles[1, 0] - kap*push_2[0]],
                [state.spindle_poles[1, 1], state.spindle_poles[1, 1] - kap*push_2[1]],
                'tab:purple', linewidth=4, label=f'F_push2 = {LA.norm(push_2):.2f}')
        ax.plot([state.spindle_poles[1, 0], state.spindle_poles[1, 0] + kap*pull_2[0]],
                [state.spindle_poles[1, 1], state.spindle_poles[1, 1] + kap*pull_2[1]],
                'tab:blue', linewidth=4, label=f'F_pull2 = {LA.norm(pull_2):.2f}')
        
        # Net forces on each pole
        ax.plot([state.spindle_poles[0, 0], state.spindle_poles[0, 0] + kap*force1[0]],
                [state.spindle_poles[0, 1], state.spindle_poles[0, 1] + kap*force1[1]],
                'tab:orange', ls='--', linewidth=4, 
                label=f'F_pole1 = {LA.norm(force1):.3f} pN', zorder=90)
        ax.plot([state.spindle_poles[1, 0], state.spindle_poles[1, 0] + kap*force2[0]],
                [state.spindle_poles[1, 1], state.spindle_poles[1, 1] + kap*force2[1]],
                'k', ls='--', linewidth=4, 
                label=f'F_pole2 = {LA.norm(force2):.3f} pN', zorder=90)
        
        # Net force
        com = np.mean(state.spindle_poles, axis=0)
        net_force = force1 + force2
        if LA.norm(net_force) > 0:
            f_net_unit = net_force / LA.norm(net_force)
        else:
            f_net_unit = np.zeros(2)
        
        ax.plot([com[0], com[0] + kap*f_net_unit[0]],
                [com[1], com[1] + kap*f_net_unit[1]],
                'k', linewidth=4, label=f'Total force = {LA.norm(net_force):.2f} pN')
    
    # Text annotations - positioned based on cell type
    mean_length = np.mean(state.mt_lengths)
    
    plt.text(right - 0.75, top - 0.15, 
             f'Total force = {LA.norm(force1 + force2):.3f}', fontsize=10)
    plt.text(right - 0.75, top - 0.25, 
             f'F_pull / F_push = {ratio:.3f}', fontsize=10)
    plt.text(right - 0.75, top - 0.35, 
             f'N_pull = {np.sum(count_pulling)}, N_push = {np.sum(count_pushing)}', 
             fontsize=10)
    
    # Color indicator for pull/push ratio
    rect = patches.Rectangle((right - 0.6, top - 0.85), 0.15, 0.15, 
                             facecolor=create_color_gradient(ratio/100))
    ax.add_patch(rect)
    
    # Chromosomes - different handling for different cell types
    com = np.mean(state.spindle_poles, axis=0)
    
    if isinstance(config, CElegansPNCConfig):
        # For PNC, show circular envelope instead of chromosomes
        theta = np.linspace(0, 2*np.pi, 360)[:-1]
        
        spindle_envelope = np.zeros((359, 2))
        spindle_envelope[:, 0] = config.spindle_half_length * np.cos(theta)
        spindle_envelope[:, 1] = config.spindle_half_length * np.sin(theta)
        
        # Rotate around origin
        R = np.array([
            [np.cos(state.spindle_angle), -np.sin(state.spindle_angle)],
            [np.sin(state.spindle_angle),  np.cos(state.spindle_angle)]
        ])
        spindle_envelope = spindle_envelope @ R.T
        
        # Translate to COM
        spindle_envelope = spindle_envelope + com
        
        # Plot envelope
        ax.plot(spindle_envelope[:, 0], spindle_envelope[:, 1], 
                linewidth=7, color='grey', zorder=1)
        plt.fill(spindle_envelope[:, 0], spindle_envelope[:, 1], 
                color='lightsteelblue', edgecolor='darkgrey', linewidth=2)
        
        # Add line connecting poles
        ax.plot([state.spindle_poles[0, 0], state.spindle_poles[1, 0]],
                [state.spindle_poles[0, 1], state.spindle_poles[1, 1]],
                color='grey', linewidth=7, zorder=5)
    else:
        # For other cell types, draw chromosomes
        chrom_angle = np.pi/2 + state.spindle_angle
        n_chrom = 4
        w = config.spindle_width
        dd = 0.8 * w
        
        xs = np.linspace(com[0] - dd*np.cos(chrom_angle), 
                        com[0] + dd*np.cos(chrom_angle), n_chrom)
        ys = np.linspace(com[1] - dd*np.sin(chrom_angle), 
                        com[1] + dd*np.sin(chrom_angle), n_chrom)
        
        for i in range(n_chrom):
            line1, line2 = draw_chromosome(xs[i], ys[i], chrom_angle)
            plt.plot([state.spindle_poles[0, 0], xs[i]], 
                    [state.spindle_poles[0, 1], ys[i]], color='g', linewidth=10)
            plt.plot([state.spindle_poles[1, 0], xs[i]], 
                    [state.spindle_poles[1, 1], ys[i]], color='g', linewidth=10)
            plt.plot(line1[0], line1[1], color='dodgerblue', linewidth=14)
            plt.plot(line2[0], line2[1], color='dodgerblue', linewidth=14)
    
    ax.legend(ncol=1, loc='upper left')
    
    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    if cell_type=='zebrafish_endo':
        plt.savefig(output_dir / f'{config.cell_type}_{config.endo_index}_{step}.pdf', 
                    bbox_inches='tight', pad_inches=1, dpi=300)

    else:
        plt.savefig(output_dir / f'{config.cell_type}_{step}.pdf', 
                    bbox_inches='tight', pad_inches=1)
    plt.close(fig)


def create_color_gradient(value: float) -> tuple:
    """Create color gradient from green (0) to white (0.5) to red (1)."""
    value = np.clip(value, 0, 1)
    
    if value <= 0.5:
        cmap = mcolors.LinearSegmentedColormap.from_list('my_cmap', ['green', 'white'])
        norm = mcolors.Normalize(vmin=0, vmax=0.5)
    else:
        cmap = mcolors.LinearSegmentedColormap.from_list('my_cmap', ['white', 'red'])
        norm = mcolors.Normalize(vmin=0.5, vmax=1)
    
    return cmap(norm(value))


def draw_chromosome(x: float, y: float, angle: float) -> Tuple[np.ndarray, np.ndarray]:
    """Draw X-shaped chromosome."""
    rc = 0.05
    
    line1 = np.array([
        [x - rc*np.cos(angle - np.pi/4), x + rc*np.cos(angle - np.pi/4)],
        [y - rc*np.sin(angle - np.pi/4), y + rc*np.sin(angle - np.pi/4)]
    ])
    line2 = np.array([
        [x - rc*np.cos(angle + np.pi/4), x + rc*np.cos(angle + np.pi/4)],
        [y - rc*np.sin(angle + np.pi/4), y + rc*np.sin(angle + np.pi/4)]
    ])
    
    return line1, line2


# ============================================================================
# EXCEL OUTPUT
# ============================================================================

def save_results_to_excel(
    config: CellConfig,
    trajectory: List[SimulationState],
    angles: np.ndarray,
    centers: np.ndarray,
    push_forces_list: List[np.ndarray],
    pull_forces_list: List[np.ndarray],
    output_dir: Path,
    run_id: int = 1
):
    """
    Save simulation results to Excel file matching original format.
    
    Creates Excel file with separate sheets for each variable:
    - Angle, Center, Y-pos, N_pull, N_push, Ratio, Pull_t, Push_t
    Each sheet has column "Run X" where X is the run_id.
    """
    
    # Extract data from trajectory
    df_angle = angles.tolist()
    df_center = [LA.norm(center) for center in centers]  # Distance from origin
    df_ypos = centers[:, 1].tolist()  # Y-position
    
    # Count pulling and pushing MTs at each timestep
    df_count_pull = [np.sum(state.which_bind) for state in trajectory]
    df_count_push = [np.sum(state.which_push) for state in trajectory]
    
    # Calculate ratio and total forces at each timestep
    df_ratio = []
    df_pull_t = []
    df_push_t = []
    
    for i in range(len(trajectory)):
        push_force = push_forces_list[i]
        pull_force = pull_forces_list[i]
        
        pull_total = np.sum(pull_force[0] + pull_force[1], axis=0)
        push_total = np.sum(push_force[0] + push_force[1], axis=0)
        
        pull_mag = LA.norm(pull_total)
        push_mag = LA.norm(push_total)
        
        if (pull_mag + push_mag) == 0:
            ratio = 0
        else:
            ratio = 100 * pull_mag / (pull_mag + push_mag)
        
        df_ratio.append(ratio)
        df_pull_t.append(pull_mag)
        df_push_t.append(push_mag)
    
    # Create dictionary matching original format
    data_dict = {
        "Angle": df_angle,
        "Center": df_center,
        "Y-pos": df_ypos,
        "N_pull": df_count_pull,
        "N_push": df_count_push,
        "Ratio": df_ratio,
        "Pull_t": df_pull_t,
        "Push_t": df_push_t
    }
    
    # Save to Excel with original format
    excel_path = output_dir / f"{config.cell_type}_data.xlsx"
    
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        for sheet_name, data in data_dict.items():
            df = pd.DataFrame({f"Run {run_id}": data})
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    
    print(f"  - Excel: {excel_path.name}")
    return excel_path


def simulate(config: CellConfig, run_id: int = 1) -> Dict[str, Any]:
    """Run complete simulation."""
    
    print(f"\n{'='*80}")
    print(f"Starting {config.cell_type} simulation")
    print(f"{'='*80}\n")
    
    # Initialize
    np.random.seed(config.seed)
    random.seed(config.seed)
    
    cell = make_cell_boundary(config, t_time=0)
    spindle_poles, spindle_angle = initialize_spindle(config, cell)
    fgs = make_fgs(config, cell, frame=1)
    
    (astral_mts, astral_angles, mt_states, which_push, which_bind,
     free_fgs, astral_which_fg, mt_lengths) = initialize_microtubules(
        config, cell, spindle_poles, spindle_angle, fgs
    )
    
    # Initial state
    state = SimulationState(
        time=0.0,
        cell_boundary=cell,
        fgs=fgs,
        spindle_poles=spindle_poles,
        spindle_angle=spindle_angle,
        astral_mts=astral_mts,
        astral_angles=astral_angles,
        mt_states=mt_states,
        mt_lengths=mt_lengths,
        which_bind=which_bind,
        which_push=which_push,
        free_fgs=free_fgs,
        astral_which_fg=astral_which_fg
    )
    
    print(f"Initialized:")
    print(f"  Cell boundary: {len(cell)} points")
    print(f"  Force generators: {len(fgs)}")
    print(f"  Astral MTs: {config.n_astral_mts} per pole")
    print(f"  Total time: {config.total_time} s")
    print(f"  Time step: {config.time_step} s")
    print(f"  Total steps: {int(config.total_time / config.time_step)}\n")
    
    # Create output directory
    if isinstance(config, ZebrafishEndoConfig):
        output_dir = config.output_dir / f"SM_{config.cell_type}_cell_{config.endo_index}_run_{run_id}"
    else:
        output_dir = config.output_dir / f"SM_{config.cell_type}_run_{run_id}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save initial configuration
    save_config(config, output_dir / "config.json")
    
    # Storage
    trajectory = [state]
    angles = [np.rad2deg(spindle_angle)]
    centers = [np.mean(spindle_poles, axis=0)]
    push_forces_list = []
    pull_forces_list = []
    
    # Plot initial state
    push_forces, pull_forces = calculate_forces(config, state)
    push_forces_list.append(push_forces)
    pull_forces_list.append(pull_forces)
    plot_cell(config, state, push_forces, pull_forces, 0, output_dir)
    
    # Time loop
    n_steps = int(config.total_time / config.time_step)
    
    # For zebrafish: determine cell update interval
    if isinstance(config, ZebrafishEndoConfig):
        frame_rate = config.extras["frame_rate"]
        update_interval = int(frame_rate / config.time_step)
    else:
        update_interval = None
    
    for step in range(n_steps):
        t_time = (step + 1) * config.time_step
        
        # Update cell shape for zebrafish at frame intervals
        if isinstance(config, ZebrafishEndoConfig) and update_interval is not None:
            if (step + 1) % update_interval == 0 and t_time < config.total_time:
                print(f"  Updating cell shape at t={t_time:.2f}s (step {step+1})")
                
                # Update cell boundary
                new_cell = update_cell_zebrafish(config, t_time)
                
                # Update FGs on new cell boundary
                new_fgs = make_fgs(config, new_cell, state.spindle_poles, 
                                   frame=int(t_time / frame_rate) + 1)
                
                # Transfer binding information to new FG array
                new_free_fgs = transfer_free_fgs(len(new_fgs), state.free_fgs)
                new_astral_which_fg = transfer_astral_which_fg(len(new_fgs), state.astral_which_fg)
                
                # Update state with new geometry
                state = SimulationState(
                    time=state.time,
                    cell_boundary=new_cell,
                    fgs=new_fgs,
                    spindle_poles=state.spindle_poles,
                    spindle_angle=state.spindle_angle,
                    astral_mts=state.astral_mts,
                    astral_angles=state.astral_angles,
                    mt_states=state.mt_states,
                    mt_lengths=state.mt_lengths,
                    which_bind=state.which_bind,
                    which_push=state.which_push,
                    free_fgs=new_free_fgs,
                    astral_which_fg=new_astral_which_fg,
                    velocity_com=state.velocity_com
                )
        
        # Calculate forces
        push_forces, pull_forces = calculate_forces(config, state)
        
        # Update spindle
        new_poles, new_angle, velocity = update_spindle_position(
            config, state, push_forces, pull_forces
        )
        
        delta_angle = new_angle - state.spindle_angle
        
        # Update state with new spindle configuration
        state = SimulationState(
            time=state.time,
            cell_boundary=state.cell_boundary,
            fgs=state.fgs,
            spindle_poles=new_poles,
            spindle_angle=new_angle,
            astral_mts=state.astral_mts,
            astral_angles=state.astral_angles,
            mt_states=state.mt_states,
            mt_lengths=state.mt_lengths,
            which_bind=state.which_bind,
            which_push=state.which_push,
            free_fgs=state.free_fgs,
            astral_which_fg=state.astral_which_fg,
            velocity_com=velocity
        )
        
        # Update microtubules
        state = update_microtubules(config, state, delta_angle)
        
        # Store
        trajectory.append(state)
        angles.append(np.rad2deg(state.spindle_angle))
        centers.append(np.mean(state.spindle_poles, axis=0))
        
        # Calculate and store forces for this timestep
        push_forces, pull_forces = calculate_forces(config, state)
        push_forces_list.append(push_forces)
        pull_forces_list.append(pull_forces)
        
        # Plot at specified intervals
        if (step + 1) % config.plot_interval == 0:
            push_forces, pull_forces = calculate_forces(config, state)
            plot_cell(config, state, push_forces, pull_forces, step + 1, output_dir)
        
        # Progress
        if (step + 1) % 100 == 0:
            print(f"  Step {step+1}/{n_steps} - " +
                  f"t={state.time:.2f}s - " +
                  f"angle={np.rad2deg(state.spindle_angle):.1f}°")
    
    print(f"\n{'='*80}")
    print("Simulation complete!")
    print(f"{'='*80}\n")
    
    # Convert lists to arrays
    angles_array = np.array(angles)
    centers_array = np.array(centers)
    
    # Save to Excel (matching original format)
    save_results_to_excel(
        config, trajectory, angles_array, centers_array,
        push_forces_list, pull_forces_list, output_dir, run_id
    )
    
    print(f"Results saved to: {output_dir}")
    print(f"  - Plots: {len(list(output_dir.glob('*.pdf')))} frames")
    print(f"  - Data: Excel file with all results")
    
    results = {
        "config": config,
        "trajectory": trajectory,
        "angles": angles_array,
        "centers": centers_array,
        "final_state": state,
        "push_forces": push_forces_list,
        "pull_forces": pull_forces_list
    }
    
    return results


# ============================================================================
# CONFIGURATION I/O
# ============================================================================

def save_config(config: CellConfig, filepath: Path):
    """Save configuration to JSON file."""
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    config_dict = {}
    for key, value in config.__dict__.items():
        if isinstance(value, Path):
            config_dict[key] = str(value)
        elif isinstance(value, np.ndarray):
            config_dict[key] = value.tolist()
        else:
            config_dict[key] = value
    
    with open(filepath, 'w') as f:
        json.dump(config_dict, f, indent=2)
    
    print(f"Configuration saved to {filepath}")


def load_config(filepath: Path) -> CellConfig:
    """Load configuration from JSON file."""
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"Config file not found: {filepath}")
    
    with open(filepath) as f:
        data = json.load(f)
    
    # Convert path strings back to Path objects
    for key in ['data_dir', 'output_dir', 'movie_info_file', 
                'cell_images_dir', 'spindle_images_dir', 'junction_data_dir']:
        if key in data and data[key] is not None:
            data[key] = Path(data[key])
    
    # Determine cell type and create appropriate config
    cell_type = data.get('cell_type')
    
    if cell_type == 'follicle_epithelial':
        return FollicleEpithelialConfig(**data)
    elif cell_type == 'neuroblast':
        return NeuroblastConfig(**data)
    elif cell_type == 'celegans_pnc':
        return CElegansPNCConfig(**data)
    elif cell_type == 'celegans_spindle':
        return CElegansSpindleConfig(**data)
    elif cell_type == 'zebrafish_endo':
        return ZebrafishEndoConfig(**data)
    else:
        raise ValueError(f"Unknown cell_type: {cell_type}")


# ============================================================================
# EXAMPLE CONFIGURATIONS
# ============================================================================

def example_follicle():
    """Example: Follicular epithelial cell."""
    config = FollicleEpithelialConfig(
        time_step=0.05,
        total_time=30.0,
        spindle_half_length=0.4,
        spindle_width=0.32,
        spindle_angle_init=np.deg2rad(90),
        n_astral_mts=100,
        fg_density=100,
        mt_mean_length=1.0,
        mt_length_stdev=1.0,
        pull_force=5.0,
        push_force=0.0
    )
    
    results = simulate(config, run_id=1)
    return results


def example_neuroblast():
    """Example: Neuroblast with apical FGs."""
    config = NeuroblastConfig(
        time_step=0.05,
        total_time=30.0,
        spindle_half_length=0.4,
        spindle_width=0.32,
        spindle_angle_init=np.deg2rad(90),
        n_astral_mts=100,
        fg_density_basal=100,
        fg_density_apical=100,
        mt_mean_length=1.0,
        mt_length_stdev=1.0,
        pull_force=5.0,
        push_force=0.0
    )
    
    results = simulate(config, run_id=1)
    return results


def example_celegans_pnc():
    """Example: C. elegans PNC."""
    config = CElegansPNCConfig(
        time_step=0.05,
        total_time=30.0,
        spindle_half_length=0.5,
        spindle_width=0.5,
        spindle_angle_init=0,  # Auto-set to 90°
        n_astral_mts=100,
        fg_density_anterior=40,
        fg_density_posterior=60,
        mt_mean_length=9.0,
        mt_length_stdev=0.167,
        pull_force=5.0,
        push_force=0.0
    )
    
    results = simulate(config, run_id=1)
    return results


def example_celegans_spindle():
    """Example: C. elegans spindle positioning."""
    config = CElegansSpindleConfig(
        time_step=0.05,
        total_time=10.0,
        spindle_half_length=0.9,
        spindle_width=0.5,
        spindle_angle_init=0,
        n_astral_mts=100,
        fg_density_anterior=40,
        fg_density_posterior=60,
        mt_mean_length=9.0,
        mt_length_stdev=0.167,
        pull_force=5.0,
        push_force=5.0
    )
    
    results = simulate(config, run_id=1)
    return results


def example_zebrafish_uniform():
    """Example: Zebrafish with uniform FG distribution."""
    config = ZebrafishEndoConfig(
        time_step=0.05,
        total_time=30.0,# If commented out, it will run full simulation
        spindle_half_length=0.18,
        spindle_width=0.08,
        spindle_angle_init=0,
        n_astral_mts=100,
        fg_density=10,
        mt_mean_length=1.0,
        mt_length_stdev=1.0,
        pull_force=5.0,
        push_force=0.0,
        endo_index=11,
        fg_distribution="uniform",
        movie_info_file=Path("./data/Movie_info.xlsx"),
        cell_images_dir=Path("./data/cells"),
        spindle_images_dir=Path("./data/spindles")
    )
    
    results = simulate(config, run_id=1)
    return results


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    import sys
    
    print("=" * 80)
    print("Spindle Positioning Simulation - Refactored Version")
    print("=" * 80)
    print()
    
    # Get run_id from SLURM array task ID if available, otherwise default to 1
    run_id = int(os.getenv('SLURM_ARRAY_TASK_ID', 1))
    
    # Check command line arguments
    if len(sys.argv) > 1:
        config_file = Path(sys.argv[1])
        if config_file.exists():
            print(f"Loading configuration from: {config_file}")
            config = load_config(config_file)
            results = simulate(config, run_id)
            
            # Save results
            output_dir = config.output_dir / f"run_{config.cell_type}"
            print(f"\nResults saved to: {output_dir}")
        else:
            print(f"Error: Config file not found: {config_file}")
            sys.exit(1)
    else:
        print("No config file provided. Running default example (C. elegans spindle).")
        print()
        results = example_zebrafish_uniform()
    
    print()
    print("=" * 80)
    print("Simulation complete!")
    print("=" * 80)

#!/usr/bin/env python
# coding: utf-8
"""
Unified spindle-orientation simulation model.

Combines three previously separate simulation scripts, developed for:
  - Follicular epithelial cells and neuroblasts (FE / NB)
  - Zebrafish embryonic cells (spindle orientation; 'junc' and 'ever'/UF FG modes)
  - C. elegans embryos (PNC: pronuclei centering, and metaphase/anaphase spindle)

into a single file for public release. Functions that are conceptually the
same across cell types but were implemented separately are kept as separate,
clearly-named functions (e.g. make_fgs_follicle / make_fgs_celegans_pnc /
make_fgs_zebrafish_junc etc.). Functions that were identical or nearly
identical across all three source scripts are kept as single shared
functions, with small cell-type/mode branches inlined where they diverged
only in a few places.
"""

import numpy as np
import matplotlib.pyplot as plt
import cv2
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import pandas as pd
import openpyxl
import os, re
from datetime import date
import time, sys
import math
import io
import contextlib

# Directory holding the zebrafish/"endo" path's external tracking data:
# cells_data.xlsx, Movie_info.xlsx, cells/cell_<c>/..., and
# spindles/spindle_<c>/... . Ships as a "data" folder in the repo root.
DATA_DIR = 'data'

# Toggles for expensive per-run I/O (per-frame PDF plots, and the raw per-MT
# xlsx log). Both default True so a plain `python spindle_model.py
# ...` (Slurm) invocation keeps producing every artifact exactly as before.
# run_simulation() (see the notebook/UI section near the end of this file)
# sets these False by default for fast interactive runs.
SAVE_FRAME_PLOTS = True
SAVE_RAW_MT_LOG = True
import pickle
import random
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from numpy import linalg as LA
from scipy.stats import norm
import shapely
import shapely.geometry
from shapely import box, LineString, normalize, Polygon
from shapely.geometry import LineString, Polygon, Point
from shapely.geometry.polygon import orient
from shapely.ops import unary_union
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
from matplotlib.path import Path
from scipy.optimize import root


# ============================================================
# SHARED / CELL-TYPE-AGNOSTIC FUNCTIONS
# (identical, or effectively identical, across all three original scripts)
# ============================================================
def angle_between_vectors(v1, v2):
    dot_product = np.dot(v1, v2)
    magnitude_v1 = np.linalg.norm(v1)
    magnitude_v2 = np.linalg.norm(v2)
    cos_angle = dot_product / (magnitude_v1 * magnitude_v2)
    angle = np.arccos(cos_angle)
    angle_degrees = np.degrees(angle)
    return angle_degrees


# def slipping_f(i,j,a,b,cell, astral_MTs,pole,orig_length,spindle_angle,n_astro,cut_part):
#     old_astral_MTs=astral_MTs
#     #print(f'orig length={orig_length}, norm = {LA.norm(pole-astral_MTs[-1])}')
#     astral_angles=np.zeros((2,n_astro))
#     astral_angles[0]=np.linspace(spindle_angle-spread/2,spindle_angle+spread/2, n_astro)
#     astral_angles[1]=np.linspace(spindle_angle+np.pi-spread/2,spindle_angle+np.pi+spread/2, n_astro)
# #     og_angle=astral_angles[i,j]
#     dist=distance_matrix(np.array([astral_MTs[-1]]),cell)[0]
#     k = 2
#     result = np.argpartition(dist, k)
#     cell_small1,cell_small2=result[:k]
#     small_dist1,small_dist2=dist[result[:k]]
#     #Find tangent vector
#     c=cell[cell_small1]-cell[cell_small1-1] #Tangent is always points counterclockwise from i to i+1
#     c=c/LA.norm(c)
#     #Find the normal vector
#     n=-normal_vector_to_ellipse(a,b, astral_MTs[-1])
#     #Find astral vector
#     a_vec=astral_MTs[-1]-astral_MTs[0] #astral vector
#     F_push=min(push, config.EI*np.pi*np.pi/LA.norm(a_vec)/LA.norm(a_vec))
#     angle = math.atan2(np.linalg.det([n,a_vec]),np.dot(n,a_vec))
#     slip=(F_push/config.mu_fric)*abs(np.sin(angle))*time_step + cut_part#total slip distance

#     #slip=max(growth_rate*time_step,cut_part) #Slip happens due to MT growth and subsequent pushing, it cant be greater than max growth amount. Exception: Pole moves towards cortex
#     # print(f'slip={slip}')
#     while (slip>=0.01*config.natural_spacing):
#         # print(f'slip I={i}')
#         dist=distance_matrix(np.array([astral_MTs[-1]]),cell)[0]
#         k = 2
#         result = np.argpartition(dist, k)
#         cell_small1,cell_small2=result[:k]
#         small_dist1,small_dist2=dist[result[:k]]
#         #Find tangent vector
#         c=cell[cell_small1]-cell[cell_small1-1] #Tangent is always points counterclockwise from i to i+1
#         c=c/LA.norm(c)
#         #Find the normal vector
#         n=-normal_vector_to_ellipse(a,b, astral_MTs[-1])
#         #Find astral vector
#         a_vec=astral_MTs[-1]-astral_MTs[0] #astral vector
#         """
#         If angle between a_vec and c is acute, it means slip is along c. If not, the slip is along -c.
#         """
#         #determine angle between astral vector and tangent
#         ac_angle=angle_between_vectors(a_vec,c)
#         if (ac_angle>90):
#             slip_t_unit=-c
#         else:
#             slip_t_unit=c

#         astral_MTs[-1]=astral_MTs[-1]+config.natural_spacing*slip_t_unit
#         dummy=astral_MTs
#         # print("after slip", LA.norm(astral_MTs[-1]-pole))
#         virtual, _ =intersect_cell(a,b,get_astral_angle(dummy), pole, cell)
#         if (LA.norm(virtual-pole)>orig_length+config.natural_spacing):
#             dist=distance_matrix(np.array([astral_MTs[-1]]),cell)[0]
#             k = 2
#             result = np.argpartition(dist, k)
#             cell_small1,cell_small2=result[:k]
#             small_dist1,small_dist2=dist[result[:k]]
#             # Closest point on a cell
#             astral_MTs[-1]=cell[cell_small1]

#             # print("after tang", LA.norm(astral_MTs[-1]-pole))
#         else:
#             astral_MTs[-1], _ =intersect_cell(a,b,get_astral_angle(dummy), pole, cell)
#             # print("after intersect", LA.norm(astral_MTs[-1]-pole))
#         slip=slip-1*config.natural_spacing
#         i=i+1
#         if (i>20):
#             break
#         # print(f'SLIP FINAL={LA.norm(astral_MTs[-1]-pole)}, coords={astral_MTs[-1]}')
#     return astral_MTs


def bounded_normal_random(mean, stdev):
    while True:
        num = random.gauss(mean, stdev)  # Generate from normal distribution
        if 0 <= num <= 1:  # Check if within bounds
            return num


def calculate_partial_perimeter(coordinates, i, j):
    # Ensure the input is a NumPy array
    coordinates = np.asarray(coordinates)

    # Check that the array is not empty and has the correct shape
    if coordinates.size == 0 or coordinates.shape[1] != 2:
        raise ValueError("Input must be a non-empty array with shape (N, 2)")

    # Normalize indices to ensure i < j
    if i > j:
        i, j = j, i

    # Extract the subarray of coordinates between indices i and j (inclusive)
    subarray = coordinates[i:j + 1]

    # Calculate the distance between consecutive points
    distances = np.sqrt(np.sum(np.diff(subarray, axis=0) ** 2, axis=1))

    # Add the distance between the last and the first point to close the loop
    closing_distance = np.sqrt(np.sum((subarray[0] - subarray[-1]) ** 2))

    # Sum all distances to get the perimeter
    perimeter = np.sum(distances) + closing_distance

    return perimeter


def calculate_perimeter(coordinates):
    # Ensure the input is a NumPy array
    coordinates = np.asarray(coordinates)

    # Check that the array is not empty and has the correct shape
    if coordinates.size == 0 or coordinates.shape[1] != 2:
        raise ValueError("Input must be a non-empty array with shape (N, 2)")

    # Calculate the distance between consecutive points
    distances = np.sqrt(np.sum(np.diff(coordinates, axis=0) ** 2, axis=1))

    # Add the distance between the last and the first point to close the loop
    closing_distance = np.sqrt(np.sum((coordinates[0] - coordinates[-1]) ** 2))

    # Sum all distances to get the perimeter
    perimeter = np.sum(distances) + closing_distance

    return perimeter


def check_bind(i, j, astral, spindle_poles, spots, free_spots, astral_which_spot):
    """
    Find the distance between the astral MT tip and the closest FG. If they are close enough, they bind together.
    """
    dist = distance_matrix(np.array([astral]), spots)
    if (np.min(dist[0]) <= config.max_interact_dist and free_spots[np.argmin(dist[0])] == 0 and random.uniform(0,
                                                                                                               1) <= prob_dyn_bind):
        free_spots[np.argmin(dist[0])] = 1
        bind = 1
        astral_which_spot[np.argmin(dist[0]), 0] = i
        astral_which_spot[np.argmin(dist[0]), 1] = j
    else:
        bind = 0
    return bind, free_spots, astral_which_spot


def check_push(a, b, astral, state, astral_angles, spindle_poles, cell):
    """
    Determines if a given MT is long enough to push
    or not.
    """
    end, _ = intersect_cell(a, b, astral_angles, spindle_poles, cell)

    # if (bind==0 and abs(LA.norm(end)-LA.norm(astral))<=time_step*growth_rate/50 and state==1):

    # math.isclose(abs(LA.norm(end)-LA.norm(astral)), 0, abs_tol=1e-12)
    # if (bind==0 and abs(LA.norm(end)-LA.norm(astral))<=min_push_dist and state==1):
    if (math.isclose(abs(LA.norm(end) - LA.norm(astral)), 0, abs_tol=1e-15) and state == 1):
        push = 1
    else:
        push = 0

    return push


def check_push_bind_init(i, j, astral, spindle_poles, spots, free_spots, astral_which_spot, state, astral_angles, cell):
    """
    Combined function that checks both binding and pushing conditions for an astral MT. Probability to bind is 1.
    Reordered workflow: first checks proximity to cell end, then binding distance.

    Parameters:
        i, j: MT indices
        astral: Astral MT tip coordinates
        spindle_poles: Spindle poles coordinates
        spots: FG spot coordinates
        free_spots: Array indicating available spots
        astral_which_spot: Tracking array
        state: MT state (1=growth, -1=shrink)
        astral_angles: MT angles
        cell: Cell geometry

    Returns:
        bind: 1 if bound, 0 otherwise
        push: 1 if pushing, 0 otherwise
        free_spots: Updated free spots array
        astral_which_spot: Updated tracking array
    """
    # Initialize outputs
    bind = 0
    push = 0

    dist = distance_matrix(np.array([astral]), spots)
    # check if the MT is close enough to a FG spot
    # and if the spot is free
    if (np.min(dist[0]) <= config.max_interact_dist and free_spots[np.argmin(dist[0])] == 0):

        free_spots[np.argmin(dist[0])] = 1
        bind = 1
        state = -1
        astral_which_spot[np.argmin(dist[0]), 0] = i
        astral_which_spot[np.argmin(dist[0]), 1] = j
    else:
        end, _ = intersect_cell(a, b, astral_angles, spindle_poles[i], cell)
        # Check if the astral MT is close enough to the cortex
        if (math.isclose(abs(LA.norm(end) - LA.norm(astral)), 0, abs_tol=config.push_dist)):
            # If not binding, check if the MT is pushing
            push = 1
            state = 1

    return bind, push, state, free_spots, astral_which_spot


def check_push_junc(a, b, astral, bind, state, astral_angles, spindle_poles, cell, spots):
    """
    For endo cells, no pushing in the region of cell-cell junctions.
    """
    end, _ = intersect_cell(a, b, astral_angles, spindle_poles, cell)
    astral_angles = normalize_angles(astral_angles)

    opening_spot = spots[0] - spindle_poles
    closing_spot = spots[-1] - spindle_poles

    opening_angle = normalize_angles(np.arctan2(opening_spot[1], opening_spot[0]))

    closing_angle = normalize_angles(np.arctan2(closing_spot[1], closing_spot[0]))

    if (bind == 0 and abs(LA.norm(end) - LA.norm(astral)) <= time_step * config.growth_rate / 2 and state == 1 and (
    not opening_angle < astral_angles < closing_angle)):  # somnitel'no
        push = 1
    else:
        push = 0

    return push


def check_spindle(spindle_poles, spindle_angle, cell, r, w):
    spindle = generate_spindle(spindle_poles, spindle_angle, r, w)
    cell = Polygon(cell)

    # dist_pole1 = distance_to_boundary(spindle_poles[0], cell)
    # dist_pole2 = distance_to_boundary(spindle_poles[1], cell)
    # min_dist = min(dist_pole1, dist_pole2)

    # We assume that spindle is an ellipse in its rigid state, but estimate drag as a for a sphere
    spindle_body = generate_spindle(spindle_poles, spindle_angle, r, w)
    spindle = Polygon(spindle_body)
    min_distance = config.min_cortex_dist
    """
    More robust version with:
    - Polygon validation
    - Edge case handling
    - Distance verification
    """
    # print(f'check spindle')
    # try:
    # cell = make_valid(Polygon(cell))
    # spindle = make_valid(Polygon(spindle_body))
    # print(f'check spindle TRY')
    # Verify minimum distance
    if cell.boundary.distance(spindle) < min_distance:
        # print(f"Debug: Boundary condition violated! - Spindle too close to cell boundary: {cell.boundary.distance(spindle)} < {min_distance}")
        return False

    # Verify complete containment
    # print(f'cell.buffer(-min_distance).contains(spindle)={cell.buffer(-min_distance).contains(spindle)}')
    return cell.buffer(-min_distance).contains(spindle)


def close_contours(contours):
    """
    Turn a contour array into a closed loop by adding points between the start and end
    such that the distance between points is less than or equal to `min_spacing`.

    Parameters:
    -----------
    contours : np.ndarray
        The array of contour coordinates with shape (N, 2).
    min_spacing : float
        The minimal spacing between the points to be added.

    Returns:
    --------
    closed_contours : np.ndarray
        The array of contour coordinates with the newly added points to close the loop.
    """
    # Extract the start and end points
    end_point = contours[0]
    start_point = contours[-1]

    # Calculate the Euclidean distance between the start and end points
    distance = np.linalg.norm(end_point - start_point)

    min_spacing = 0.005
    # Calculate how many points to add based on the minimal spacing
    num_points_to_add = int(np.ceil(distance / min_spacing))

    # Create the interpolation for the required number of points
    x_new = np.linspace(start_point[0], end_point[0], num_points_to_add + 2)  # +2 includes start and end points
    y_new = np.linspace(start_point[1], end_point[1], num_points_to_add + 2)

    # Remove the first and last points from the new points to avoid duplication with original points
    x_new = x_new[1:-1]
    y_new = y_new[1:-1]

    # Combine the new points into a new array
    new_points = np.column_stack((x_new, y_new))

    # Combine original contours with new points to form a closed contour
    closed_contours = np.vstack([contours, new_points])

    # Close the loop by adding the first point at the end
    # closed_contours = np.vstack([closed_contours, closed_contours[0]])

    return closed_contours


# def compute_resistance_functions(a, b):
#     """
#     Computes the translational and rotational resistance functions (Table 3.4).
#     Now handles circular cross-sections (a == b) with analytical solutions.
#
#     Parameters:
#     a : float - Semi-major axis
#     b : float - Semi-minor axis
#
#     Returns:
#     X_A, Y_A, X_C, Y_C : float - Resistance functions for translation & rotation
#     """
#
#     # Check for circular case first
#     if a == b:
#         # Analytical solutions for circle (all equal)
#         resistance = 16 / 3
#         return resistance, resistance, resistance, resistance
#
#     # Original check for invalid input
#     if a < b:
#         raise ValueError("Semi-major axis (a) must be greater than or equal to semi-minor axis (b).")
#
#     # Compute eccentricity
#     e = np.sqrt(1 - (b ** 2 / a ** 2))
#
#     # Compute logarithmic correction term
#     L = np.log((1 + e) / (1 - e))
#
#     # Compute resistance functions
#     X_A = (8 / 3) * (e ** 3 / (-2 * e + (1 + e ** 2) * L))
#     Y_A = (16 / 3) * (e ** 3 / (2 * e + (1 + e ** 2) * L))
#     X_C = (4 / 3) * (e ** 3 * (1 - e ** 2) / (2 * e - (1 - e ** 2) * L))
#     Y_C = (4 / 3) * (e ** 3 * (2 - e ** 2) / (-2 * e + (1 + e ** 2) * L))
#
#     return X_A, Y_A, X_C, Y_C


# def generate_spindle(spindle_poles, spindle_angle, r, w):
#     spindle=np.zeros((359,2))
#     theta = np.linspace(0, 2*np.pi, 360)[:-1].copy()

#     com=np.array([(spindle_poles[0,0]+spindle_poles[1,0])/2,(spindle_poles[0,1]+spindle_poles[1,1])/2]) #rotation is about com(centre of mass)
#     spindle[:,0]=r*np.cos(theta)+com[0]
#     spindle[:,1]=w*np.sin(theta)+com[1]

#     spindle = rotate_points(spindle, spindle_angle) #rotate spindle to the right angle

#     return spindle


def create_color_gradient(value):
    """
    Create a color gradient from green to white to red based on the provided value.
    The value should range from 0 to 1.
    """
    if value < 0:
        value = 0
    elif value > 1:
        value = 1

    if value <= 0.5:
        cmap = mcolors.LinearSegmentedColormap.from_list('my_cmap', ['green', 'white'])
        norm = mcolors.Normalize(vmin=0, vmax=0.5)
    else:
        cmap = mcolors.LinearSegmentedColormap.from_list('my_cmap', ['white', 'red'])
        norm = mcolors.Normalize(vmin=0.5, vmax=1)

    return cmap(norm(value))


def create_simulation_directory(test_folder_path, name, task_id):
    folder_name = f"{name}_{task_id}"
    new_dir_path = os.path.join(test_folder_path, folder_name)
    # makedirs(..., exist_ok=True) instead of mkdir: on a cluster the array
    # job always gets a fresh folder_name per task_id, but a local user
    # re-running the same name locally would otherwise hit FileExistsError.
    os.makedirs(new_dir_path, exist_ok=True)
    return folder_name, new_dir_path


def delta_theta(a, b, theta):
    if a * b == 3:  # or theta <0 or theta>2*np.pi:
        # return 'What the fuck....'
        c = 1 + 3
    else:
        if (theta < 0):
            theta = theta + 2 * np.pi
        if (theta > 2 * np.pi):
            theta = theta - 2 * np.pi
        if a >= b:
            if theta <= np.pi / 2:
                return -np.arccos((a * np.cos(theta) ** 2 + b * np.sin(theta) ** 2) / np.sqrt(
                    (a * np.cos(theta)) ** 2 + (b * np.sin(theta)) ** 2)) + theta
            elif theta <= np.pi:
                return np.arccos((a * np.cos(theta) ** 2 + b * np.sin(theta) ** 2) / np.sqrt(
                    (a * np.cos(theta)) ** 2 + (b * np.sin(theta)) ** 2)) + theta
            elif theta <= 3 * np.pi / 2:
                return -np.arccos((a * np.cos(theta) ** 2 + b * np.sin(theta) ** 2) / np.sqrt(
                    (a * np.cos(theta)) ** 2 + (b * np.sin(theta)) ** 2)) + theta
            else:
                return np.arccos((a * np.cos(theta) ** 2 + b * np.sin(theta) ** 2) / np.sqrt(
                    (a * np.cos(theta)) ** 2 + (b * np.sin(theta)) ** 2)) + theta
        else:
            if theta <= np.pi / 2:
                return np.arccos((a * np.cos(theta) ** 2 + b * np.sin(theta) ** 2) / np.sqrt(
                    (a * np.cos(theta)) ** 2 + (b * np.sin(theta)) ** 2)) + theta
            elif theta <= np.pi:
                return -np.arccos((a * np.cos(theta) ** 2 + b * np.sin(theta) ** 2) / np.sqrt(
                    (a * np.cos(theta)) ** 2 + (b * np.sin(theta)) ** 2)) + theta
            elif theta <= 3 * np.pi / 2:
                return np.arccos((a * np.cos(theta) ** 2 + b * np.sin(theta) ** 2) / np.sqrt(
                    (a * np.cos(theta)) ** 2 + (b * np.sin(theta)) ** 2)) + theta
            else:
                return -np.arccos((a * np.cos(theta) ** 2 + b * np.sin(theta) ** 2) / np.sqrt(
                    (a * np.cos(theta)) ** 2 + (b * np.sin(theta)) ** 2)) + theta


def distance_to_boundary(point, shape_coords):
    """
    Calculate the minimum distance from a point to the boundary defined by shape_coords.

    Parameters:
    point (np.array): The point (x, y) for which to calculate the distance.
    shape_coords (np.array): The coordinates of the boundary shape (N, 2).

    Returns:
    float: The minimum distance from the point to the shape's boundary.
    """
    # Calculate the Euclidean distance from the point to all points in shape_coords
    distances = np.sqrt((shape_coords[:, 0] - point[0]) ** 2 + (shape_coords[:, 1] - point[1]) ** 2)

    return np.min(distances)


def draw_vectors(ax, origin, forces, labels, colors, normalized_length=0.2, text_offset=0.3):
    """
    Draw vectors defined by their endpoints.

    Parameters:
        ax: Matplotlib axis object.
        origin: The starting point of all vectors (e.g., [x, y]).
        forces: List of force vectors, where each vector is an endpoint [x, y].
        labels: List of labels for each vector.
        colors: List of colors for each vector.
        normalized_length: Desired length of the vectors (optional).
        text_offset: Offset for text labels to avoid overlap (optional).
    """
    for i, force_end in enumerate(forces):
        # Calculate the direction vector
        dx = force_end[0] - origin[0]
        dy = force_end[1] - origin[1]

        # Normalize the vector to the desired length
        if normalized_length is not None:
            magnitude = np.sqrt(dx ** 2 + dy ** 2)
            if magnitude > 0:  # Avoid division by zero
                scale = normalized_length / magnitude
                dx *= scale
                dy *= scale

        # Calculate the endpoint of the normalized vector
        end_x = origin[0] + dx
        end_y = origin[1] + dy

        # Draw the arrow
        ax.arrow(origin[0], origin[1], dx, dy,
                 head_width=0.05, head_length=0.1, fc=colors[i], ec=colors[i])

        # Calculate text position with an offset
        text_x = end_x + dx * text_offset
        text_y = end_y + dy * text_offset

        # Add text label
        ax.text(text_x, text_y, f"{labels[i]}", fontsize=8, ha='center', color=colors[i])


def enhance_cell(cell):
    """
    The spacings between raw cell contour vertices are unequal. This function fixes it.
    """
    new_cell = []
    max_interj_dist = 0.005
    # max_interj_dist=config.natural_spacing
    for i in range(len(cell)):
        if (i == len(cell) - 1):
            next = 0
        else:
            next = i + 1
        dist = LA.norm([cell[i, 0] - cell[next, 0], cell[i, 1] - cell[next, 1]])
        if (dist > max_interj_dist):
            new_cell.append(cell[i])
            spread = np.linspace(0, dist, int(dist / max_interj_dist))[1:-1]
            vec = np.array([cell[next, 0] - cell[i, 0], cell[next, 1] - cell[i, 1]]) / dist  # +cell[i]
            for j in range(len(spread)):
                new_add = spread[j] * vec + cell[i]
                new_cell.append(new_add)
        else:
            new_cell.append(cell[i])
    new_cell_ar = np.zeros((len(new_cell), 2))
    for i in range(len(new_cell)):
        new_cell_ar[i] = new_cell[i]
    return new_cell_ar


# Log function to non-linearly distribute FGs

def find_number_after_cell(input_string):
    # Regular expression to find a number after "cell"
    match = re.search(r'endo_(\d+)', input_string)
    if match:
        number = match.group(1)
    return int(number)


def find_tangential_and_normal_vectors(spindle_poles, spindle_angle, r, w, cell):
    """
    Finds tangential and normal unit vectors at the closest point of the spindle to the cell boundary.
    - spindle_poles: (2,2)
    - spindle_angle: float
    - r, w: ellipse axes
    - cell: (N,2) polygon
    """
    spindle = generate_spindle(spindle_poles, spindle_angle, r, w)
    # Find the closest point on spindle to cell
    D = distance_matrix(spindle, cell)
    idx = np.argmin(D)
    row, col = np.unravel_index(idx, D.shape)
    close_spindle = spindle[row]
    close_cell = cell[col]

    # Find which edge on the cell this is (between col and col+1)
    N = len(cell)
    segs = np.column_stack([cell, np.roll(cell, -1, axis=0)])  # N x 4
    # Find distance from close_spindle to each segment
    min_dist = np.inf
    best_i = None
    for i in range(N):
        a, b = segs[i, :2], segs[i, 2:]
        ab = b - a
        ap = close_spindle - a
        t = np.clip(np.dot(ap, ab) / np.dot(ab, ab), 0, 1)
        closest = a + t * ab
        d = LA.norm(close_spindle - closest)
        if d < min_dist:
            min_dist = d
            best_i = i
            tangent_vec = ab / LA.norm(ab)
    tangential = tangent_vec
    normal = np.array([-tangential[1], tangential[0]])  # Outward for CCW polygons

    return tangential, normal


# def slide_and_rotate_spindle(spindle_poles, spindle_angle, V, Omega, cell, r, w, config, dt):
#     """
#     Moves and rotates the spindle as far as allowed by cortex in this time step.
#     1. Try maximal safe translation along tangential direction.
#     2. Then try maximal safe rotation at the new location.
#     Debug lines are included for tracking steps.
#     If no move is allowed, returns original poles and angle.
#     """

#     print("\n[DEBUG] Starting slide_and_rotate_spindle")
#     print(f"[DEBUG] Input spindle_angle: {np.rad2deg(spindle_angle):.2f} deg, V: {V}, Omega: {Omega}, dt: {dt}")

#     # --- Translation first ---
#     tangential_vec, normal_vec = find_tangential_and_normal_vectors(spindle_poles, spindle_angle, r, w, cell)
#     V_tan = np.dot(V, tangential_vec) * tangential_vec
#     v_tan_norm = np.linalg.norm(V_tan)
#     print(f"[DEBUG] Tangential velocity component: {V_tan}, norm: {v_tan_norm}")

#     if v_tan_norm == 0:
#         print("[DEBUG] No tangential velocity; translation skipped.")
#         spindle_poles_trans = spindle_poles.copy()
#         actual_step = 0
#     else:
#         # Cap move by max_step
#         step_size = min(v_tan_norm * dt, config.max_step)
#         print(f"[DEBUG] Capped step size: {step_size} (max allowed: {config.max_step})")

#         # Binary search for largest allowed translation
#         low, high = 0, 1
#         for j in range(10):
#             frac = (low + high) / 2
#             trial_poles = spindle_poles + frac * V_tan / v_tan_norm * step_size
#             is_valid = check_spindle(trial_poles, spindle_angle, cell, r, w)
#             if is_valid:
#                 low = frac
#             else:
#                 high = frac
#             print(f"[DEBUG] Trans search iter {j}: frac={frac:.4f}, valid={is_valid}")
#         spindle_poles_trans = spindle_poles + low * V_tan / v_tan_norm * step_size
#         actual_step = low * step_size
#         print(f"[DEBUG] Final translation fraction: {low:.4f}, actual step: {actual_step:.4f} um")

#     # --- Now, rotation at new translated position ---
#     dtheta = Omega * dt
#     print(f"[DEBUG] Intended rotation (radians): {dtheta:.4f} ({np.rad2deg(dtheta):.2f} deg)")
#     low, high = 0, 1
#     for k in range(10):
#         frac = (low + high) / 2
#         trial_angle = spindle_angle + frac * dtheta
#         com = np.mean(spindle_poles_trans, axis=0)
#         R = np.array([
#             [np.cos(trial_angle - spindle_angle), -np.sin(trial_angle - spindle_angle)],
#             [np.sin(trial_angle - spindle_angle),  np.cos(trial_angle - spindle_angle)]
#         ])
#         trial_poles_rot = (spindle_poles_trans - com) @ R.T + com
#         is_valid = check_spindle(trial_poles_rot, trial_angle, cell, r, w)
#         if is_valid:
#             low = frac
#         else:
#             high = frac
#         print(f"[DEBUG] Rot search iter {k}: frac={frac:.4f}, angle={np.rad2deg(trial_angle):.2f} deg, valid={is_valid}")
#     final_angle = spindle_angle + low * dtheta
#     actual_rot = low * dtheta
#     com = np.mean(spindle_poles_trans, axis=0)
#     R = np.array([
#         [np.cos(final_angle - spindle_angle), -np.sin(final_angle - spindle_angle)],
#         [np.sin(final_angle - spindle_angle),  np.cos(final_angle - spindle_angle)]
#     ])
#     final_poles = (spindle_poles_trans - com) @ R.T + com

#     print(f"[DEBUG] Final rotation fraction: {low:.4f}, actual rotation: {actual_rot:.4f} radians ({np.rad2deg(actual_rot):.2f} deg)")
#     print(f"[DEBUG] Final spindle center: {np.mean(final_poles, axis=0)}")

#     # --- SAFEGUARD: Return original if nothing succeeded ---
#     if not check_spindle(final_poles, final_angle, cell, r, w):
#         print("[WARNING] No valid translation/rotation possible; returning original spindle configuration.")
#         return spindle_poles.copy(), spindle_angle, False
#     else:
#         print(f"[DEBUG] --- End slide_and_rotate_spindle ---\n")
#         return final_poles, final_angle, True  # Return success flag

# def slide_and_rotate_spindle(spindle_poles, spindle_angle, V, Omega, cell, r, w, config, dt):
#     """
#     Attempts maximal safe translation and rotation of spindle in cell.
#     Returns new poles, angle, and a boolean (True if any move occurred, else False).
#     """
#     print("\n[DEBUG] Starting slide_and_rotate_spindle")
#     print(f"[DEBUG] Input spindle_angle: {np.rad2deg(spindle_angle):.2f} deg, V: {V}, Omega: {Omega}, dt: {dt}")

#     # --- Maximal tangential translation ---
#     tangential_vec, _ = find_tangential_and_normal_vectors(spindle_poles, spindle_angle, r, w, cell)
#     V_tan = np.dot(V, tangential_vec) * tangential_vec
#     v_tan_norm = np.linalg.norm(V_tan)
#     step_size = min(v_tan_norm * dt, config.max_step) if v_tan_norm > 0 else 0

#     print(f"[DEBUG] Tangential velocity: {V_tan}, norm: {v_tan_norm}, step_size: {step_size}")

#     # Binary search for translation
#     low, high = 0, 1
#     for j in range(10):
#         frac = (low + high) / 2
#         trial_poles = spindle_poles + frac * V_tan / v_tan_norm * step_size if v_tan_norm > 0 else spindle_poles.copy()
#         is_valid = check_spindle(trial_poles, spindle_angle, cell, r, w)
#         if is_valid:
#             low = frac
#         else:
#             high = frac
#         print(f"[DEBUG] Trans iter {j}: frac={frac:.4f}, valid={is_valid}")
#     spindle_poles_trans = spindle_poles + low * V_tan / v_tan_norm * step_size if v_tan_norm > 0 else spindle_poles.copy()
#     moved_translation = np.any(np.abs(spindle_poles_trans - spindle_poles) > 1e-9)

#     # --- Maximal allowed rotation at translated position ---
#     dtheta = Omega * dt
#     low_rot, high_rot = 0, 1
#     for k in range(10):
#         frac = (low_rot + high_rot) / 2
#         trial_angle = spindle_angle + frac * dtheta
#         com = np.mean(spindle_poles_trans, axis=0)
#         R = np.array([
#             [np.cos(trial_angle - spindle_angle), -np.sin(trial_angle - spindle_angle)],
#             [np.sin(trial_angle - spindle_angle),  np.cos(trial_angle - spindle_angle)]
#         ])
#         trial_poles_rot = (spindle_poles_trans - com) @ R.T + com
#         is_valid = check_spindle(trial_poles_rot, trial_angle, cell, r, w)
#         if is_valid:
#             low_rot = frac
#         else:
#             high_rot = frac
#         print(f"[DEBUG] Rot iter {k}: frac={frac:.4f}, angle={np.rad2deg(trial_angle):.2f} deg, valid={is_valid}")

#     final_angle = spindle_angle + low_rot * dtheta
#     com = np.mean(spindle_poles_trans, axis=0)
#     R = np.array([
#         [np.cos(final_angle - spindle_angle), -np.sin(final_angle - spindle_angle)],
#         [np.sin(final_angle - spindle_angle),  np.cos(final_angle - spindle_angle)]
#     ])
#     final_poles = (spindle_poles_trans - com) @ R.T + com
#     moved_rotation = (abs(final_angle - spindle_angle) > 1e-9)

#     # --- Determine if any change occurred and return ---
#     has_changed = (moved_translation or moved_rotation)

#     if not check_spindle(final_poles, final_angle, cell, r, w):
#         print("[WARNING] No valid translation/rotation possible; returning original spindle configuration.")
#         return spindle_poles.copy(), spindle_angle, False
#     elif has_changed:
#         print("[DEBUG] --- End slide_and_rotate_spindle: MOVE occurred ---\n")
#         return final_poles, final_angle, True
#     else:
#         print("[DEBUG] --- End slide_and_rotate_spindle: NONE ---\n")
#         return spindle_poles.copy(), spindle_angle, False


def flip_angles(df_angle):
    if (df_angle[-1] > 90):
        for i in range(len(df_angle)):
            dummy = df_angle[i]
            df_angle[i] = 180 - dummy
    elif (df_angle[-1] < -90):
        for i in range(len(df_angle)):
            dummy = df_angle[i]
            df_angle[i] = 180 + dummy
    return df_angle


def flip_angles_NB(df_angle):
    if (df_angle[1] > 160):
        for i in range(len(df_angle)):
            dummy = df_angle[i]
            df_angle[i] = 180 - dummy
    elif (df_angle[-1] < -90):
        for i in range(len(df_angle)):
            dummy = df_angle[i]
            df_angle[i] = 180 + dummy
    return df_angle


def gaussian_spaced_indices(x, y, n_samples, mean=None, std_dev=None):
    """
    Generate indices spaced according to a Gaussian distribution between x and y.

    Args:
        x (int): Start index (left bound).
        y (int): End index (right bound).
        n_samples (int): Number of indices to generate.
        mean (float): Center of the Gaussian (default: midpoint between x and y).
        std_dev (float): Spread of the Gaussian (default: (y - x)/6 to cover 99.7% of range).

    Returns:
        np.ndarray: Array of indices (float or int).
    """
    if mean is None:
        mean = (x + y) / 2  # Midpoint as default mean
    if std_dev is None:
        std_dev = (y - x) / 6  # Default spread (3σ covers 99.7% of range)

    # Generate evenly-spaced quantiles (avoid extremes 0 and 1 to prevent infinite values)
    quantiles = np.linspace(0.001, 0.999, n_samples)

    # Map quantiles to Gaussian-distributed indices
    indices = norm.ppf(quantiles, loc=mean, scale=std_dev * 0.7)  # std_dev)

    # Clip to ensure indices stay within [x, y] (optional, but recommended)
    indices = np.clip(indices, x, y)

    return indices.astype(int)  # Convert to integer indices


def generate_spindle(spindle_poles, spindle_angle, r, w):
    theta = np.linspace(0, 2 * np.pi, 36, endpoint=False)
    com = np.mean(spindle_poles, axis=0)

    # Superellipse parameters
    a, b = r, w  # semi-major and semi-minor axes
    n = 1.4
    na = 2 / n

    # Create the superellipse shape
    x = (np.abs(np.cos(theta)) ** na) * a * np.sign(np.cos(theta))
    y = (np.abs(np.sin(theta)) ** na) * b * np.sign(np.sin(theta))
    spindle = np.vstack((x, y)).T + com  # shift to COM

    # Rotate around COM
    R = np.array([
        [np.cos(spindle_angle), -np.sin(spindle_angle)],
        [np.sin(spindle_angle), np.cos(spindle_angle)]
    ])
    spindle = (spindle - com) @ R.T + com

    return spindle


def get_astral_angle(astral):
    angle = np.arctan2(astral[-1, 1] - astral[0, 1], astral[-1, 0] - astral[0, 0])
    if (-np.pi <= angle <= 0):
        angle = angle + 2 * np.pi
    return angle


def get_contours(image_path):
    image = cv2.imread(image_path)
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    # Apply Gaussian Blur to reduce noise and detail
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)

    # Apply edge detection using Canny
    edges = cv2.Canny(blurred, 50, 150)

    # Find contours in the edge-detected image
    contours, hierarchy = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    smoothed_contours = []

    # Smoothen each contour
    for contour in contours:
        # Approximate the contour with a smoother curve
        epsilon = 0.0001 * cv2.arcLength(contour, True)
        smoothed_contour = cv2.approxPolyDP(contour, epsilon, True)
        smoothed_contours.append(smoothed_contour)

    # Concatenate all contour points into a single array
    all_points = np.concatenate([contour.squeeze() for contour in smoothed_contours])
    all_points[:, 1] = -all_points[:, 1]

    return all_points


def get_sequence(
        n_points):  # gives you the sequence of angles (theta) that acts as coordinates (x=r cos theta, y=r sin theta) where to put spots or protein machines

    sequence = np.zeros((n_points))
    step = ((np.exp(np.pi / 4)) - 1 / (np.exp(np.pi))) / n_points
    for i in range(n_points):
        sequence[i] = log_func(i * step)

    return sequence


def get_sequence_celegans(n_points):
    sequence = np.zeros((n_points) // 2 + 1)
    step = (np.exp(np.pi / 3) - 1) / (n_points // 2)
    for i in range(n_points // 2 + 1):
        sequence[i] = log_func_celegans(i * step)
    dum = -sequence[::-1]
    sequence = np.concatenate((sequence, dum[1:]))

    return sequence


def get_sequence_top(n_points):  # same but for the apical surface

    sequence = np.zeros((n_points) // 2 + 1)
    sequence2 = np.zeros((n_points) // 2)
    step = (4.67) / (n_points // 2)
    for i in range(n_points // 2):
        sequence[i] = log_func_top(i * step)
    for i in range(n_points // 2):
        sequence2[i] = log_func_top2(i * step)
    dum = sequence2[::-1]
    sequence[(n_points) // 2] = 0.5 * (sequence[-2] + dum[0])
    sequence = np.concatenate((sequence, dum))

    return sequence


def get_spindle_angle(spindle):
    angle = np.arctan2(spindle[0, 1] - spindle[1, 1], spindle[0, 0] - spindle[1, 0])
    if (-np.pi <= angle <= 0):
        angle = angle + 2 * np.pi
    return angle


def grow_astralMT(a, b, angle, C, cell, orig_length):
    # print("grow_astralMT")
    # Returns MT end point after growrth and if it touches the cortex (which_push) since state=1 because its growing
    intersect, _ = intersect_cell(a, b, angle, C, cell)

    end = C + (orig_length + config.growth_rate * time_step) * np.array([np.cos(angle), np.sin(angle)])

    """
    Checks if MT tip ends up beyond cell cortex after growth. In this case, the MT end will be set to the intersection with cell cortex
    """
    len1 = LA.norm([C[0] - intersect[0], C[1] - intersect[1]])
    len2 = LA.norm([C[0] - end[0], C[1] - end[1]])
    # if (orig_length+config.growth_rate*time_step>=config.MT_max_length):
    #     if (point_in_polygon(C+(orig_length)*np.array([np.cos(angle),np.sin(angle)]), cell)==True):
    #         return C+(orig_length)*np.array([np.cos(angle),np.sin(angle)]), 0
    #     else:
    #         return intersect,1

    if (
            len2 < len1):  # and point_in_polygon(end, cell)==True): #If grown end is closer to the pole than intersect with the cortex
        # print(f'grow end={end}, start_point={C}')
        # print(f'new MT length= {LA.norm(C-end)}')
        return end
    else:
        # print(f'grow intersect={end}, start_point={C}')
        # print(f'new MT length= {LA.norm(C-end)}')
        return intersect


def intersect_cell(a, b, angle, start_point, cell):
    """
    Deploys shapely module to find an intersection between astral MTs and cell cortex
    """
    # NumPy array
    end_point = start_point + 4 * np.array([np.cos(angle), np.sin(angle)])
    numpy_line_coords = np.array([start_point, end_point])
    # Convert to Shapely LineString
    shapely_line = shapely.geometry.LineString(numpy_line_coords)
    shapely_string = LineString(cell)
    inter = shapely.intersection(shapely_string, shapely_line)

    if (inter.geom_type == 'Point'):
        return np.array([inter.x, inter.y]), inter.geom_type
    elif (inter.geom_type == 'MultiPoint'):
        points = [p for p in inter.geoms]
        dist = []
        for i in range(len(points)):
            dist.append(LA.norm([start_point[0] - points[i].x, start_point[1] - points[i].y]))
        return np.array([points[np.argsort(dist)[0]].x, points[np.argsort(dist)[0]].y]), inter.geom_type
    else:
        # print(f'intersection type:{inter.geom_type}, start={start_point}, end={end_point}')
        # print(np.asanyarray(inter))
        intersect = intersect_cell_old(a, b, angle, start_point, cell)
        return intersect, inter.geom_type  # start_point+0.8*MT_max_length*np.array([np.cos(angle),np.sin(angle)])


def intersect_cell_old(a, b, angle, C, cell):
    """
    When astral MT's is shrinking but spindle moves towards the cortex there are 3 cases:
    1. End is outside
    2. End is inside
    This function basically tells if the end is inside keep it there, if its outside, find the cortex intersect and take the end there.
    """
    L_cell = 2 * np.pi / config.number_of_sides
    D = C
    L_min = min(a, b)  #
    j = 0
    dist = distance_matrix(np.array([C]), cell)[0]
    while (L_min > L_cell):  # and (C[0]/a)**2+(C[1]/b)**2<1.05):
        if (j > 10):
            break
        dist = distance_matrix(np.array([C]), cell)[0]
        L_min = np.min(dist)
        C = C + L_min * np.array([np.cos(angle), np.sin(angle)])
        j = j + 1

    D = C + 2 * np.array([np.cos(angle), np.sin(angle)])
    k = 2
    result = np.argpartition(dist, k)
    cell_small1, cell_small2 = result[:k]
    small_dist1, small_dist2 = dist[result[:k]]
    A = cell[cell_small1]
    B = cell[cell_small2]
    # Line AB represented as a1x + b1y = c1
    a1 = B[1] - A[1]
    b1 = A[0] - B[0]
    c1 = a1 * (A[0]) + b1 * (A[1])
    # Line CD represented as a2x + b2y = c2
    a2 = D[1] - C[1]
    b2 = C[0] - D[0]
    c2 = a2 * (C[0]) + b2 * (C[1])
    determinant = a1 * b2 - a2 * b1
    if (determinant == 0):
        return D
    else:
        x = (b2 * c1 - b1 * c2) / determinant
        y = (a1 * c2 - a2 * c1) / determinant
        intersect = np.array([x, y])
        return intersect


def log_func(x):
    return 0.5 * (np.log(x + 1 / np.exp(np.pi)))


def log_func_celegans(x):
    return -np.log(x + 1) + np.pi / 3


def log_func_top(x):
    return 0.3 * np.log(x + 1 / np.exp(1)) + np.pi / 4 + 0.3


def log_func_top2(x):
    return -0.3 * np.log(x + 1 / np.exp(1)) + 3 * np.pi / 4 - 0.3


def make_new_MT(i, j, a, b, MT, spindle_poles, astral_angles, cell, free_spots, astral_which_spot, orig_length):
    free_spots[np.argwhere((astral_which_spot[:, 0] == i) & (astral_which_spot[:, 1] == j))] = 0
    orig_length = 0
    MT[-1] = grow_astralMT(a, b, astral_angles[i, j], spindle_poles[i], cell, orig_length)

    MT = restructure(MT)

    # print(f'new MT length= {LA.norm(MT[-1]-spindle_poles[i])}')
    # print(f'new MT angle={astral_angles[i,j]}')
    return MT


def normal_vector_to_ellipse(a, b, coordinate):
    # Unpack ellipse parameters

    # Calculate the gradient of the ellipse equation
    x, y = coordinate[0], coordinate[1]
    gradient_x = (2 * x) / (a ** 2)
    gradient_y = (2 * y) / (b ** 2)

    # Normalize the gradient vector to get the unit normal vector
    magnitude = np.sqrt(gradient_x ** 2 + gradient_y ** 2)
    normal_vector = np.array([gradient_x / magnitude, gradient_y / magnitude])

    return normal_vector / LA.norm(normal_vector)


def normalize_angles(angles):
    """
    Normalize angles to the range [0, 2π].

    Parameters:
    angles: scalar, 1D array, or (2, N) array containing angles in radians

    Returns:
    normalized_angles: object of the same shape as input with angles normalized to [0, 2π]
    """
    # Convert to numpy array for ease of operations
    angles = np.asarray(angles)

    # Normalize angles to the range [0, 2π]
    normalized_angles = angles % (2 * np.pi)

    return normalized_angles


def point_in_polygon(point, polygon):
    """
    For checking if spindle poles leave the cell boundary.
    """
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


def restructure(beam):
    restructured = np.ones(np.shape(beam))
    restructured[:, 0] = np.linspace(beam[0, 0], beam[-1, 0], config.discr)
    restructured[:, 1] = np.linspace(beam[0, 1], beam[-1, 1], config.discr)
    return restructured


def rotate_and_shift_parabola(x, y, theta, xi, yi):
    # Rotation matrix
    rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                                [np.sin(theta), np.cos(theta)]])

    # Stack x and y to create coordinate pairs
    coordinates = np.vstack((x, y))

    # Apply rotation
    rotated_coordinates = rotation_matrix @ coordinates

    # Shift the coordinates by (xi, yi)
    x_shifted = rotated_coordinates[0] + xi
    y_shifted = rotated_coordinates[1] + yi

    return x_shifted, y_shifted


def sgn(x):
    """Sign function using NumPy."""
    return np.where(x > 0, 1, -1)  # Vectorized version


# def slipping_f(i, j, a, b, cell, astral_MTs, pole, orig_length, spindle_angle, n_astro, cut_part):
#     old_astral_MTs = astral_MTs
#     # print(f'orig length={orig_length}, norm = {LA.norm(pole-astral_MTs[-1])}')
#     astral_angles = np.zeros((2, n_astro))
#     astral_angles[0] = np.linspace(spindle_angle - spread / 2, spindle_angle + spread / 2, n_astro)
#     astral_angles[1] = np.linspace(spindle_angle + np.pi - spread / 2, spindle_angle + np.pi + spread / 2, n_astro)
#     #     og_angle=astral_angles[i,j]
#     dist = distance_matrix(np.array([astral_MTs[-1]]), cell)[0]
#     k = 2
#     result = np.argpartition(dist, k)
#     cell_small1, cell_small2 = result[:k]
#     small_dist1, small_dist2 = dist[result[:k]]
#     # Find tangent vector
#     c = cell[cell_small1] - cell[cell_small1 - 1]  # Tangent is always points counterclockwise from i to i+1
#     c = c / LA.norm(c)
#     # Find the normal vector
#     n = -normal_vector_to_ellipse(a, b, astral_MTs[-1])
#     # Find astral vector
#     a_vec = astral_MTs[-1] - astral_MTs[0]  # astral vector
#     F_push = min(push, config.EI * np.pi * np.pi / LA.norm(a_vec) / LA.norm(a_vec))
#     angle = math.atan2(np.linalg.det([n, a_vec]), np.dot(n, a_vec))
#     slip = (F_push / config.mu_fric) * abs(np.sin(angle)) * time_step + cut_part  # total slip distance
#
#     # slip=max(growth_rate*time_step,cut_part) #Slip happens due to MT growth and subsequent pushing, it cant be greater than max growth amount. Exception: Pole moves towards cortex
#     # print(f'slip={slip}')
#     while (slip >= 0.01 * config.natural_spacing):
#         # print(f'slip I={i}')
#         dist = distance_matrix(np.array([astral_MTs[-1]]), cell)[0]
#         k = 2
#         result = np.argpartition(dist, k)
#         cell_small1, cell_small2 = result[:k]
#         small_dist1, small_dist2 = dist[result[:k]]
#         # Find tangent vector
#         c = cell[cell_small1] - cell[cell_small1 - 1]  # Tangent is always points counterclockwise from i to i+1
#         c = c / LA.norm(c)
#         # Find the normal vector
#         n = -normal_vector_to_ellipse(a, b, astral_MTs[-1])
#         # Find astral vector
#         a_vec = astral_MTs[-1] - astral_MTs[0]  # astral vector
#         """
#         If angle between a_vec and c is acute, it means slip is along c. If not, the slip is along -c.
#         """
#         # determine angle between astral vector and tangent
#         ac_angle = angle_between_vectors(a_vec, c)
#         if (ac_angle > 90):
#             slip_t_unit = -c
#         else:
#             slip_t_unit = c
#
#         astral_MTs[-1] = astral_MTs[-1] + config.natural_spacing * slip_t_unit
#         dummy = astral_MTs
#         # print("after slip", LA.norm(astral_MTs[-1]-pole))
#         virtual, _ = intersect_cell(a, b, get_astral_angle(dummy), pole, cell)
#         if (LA.norm(virtual - pole) > orig_length + config.natural_spacing):
#             dist = distance_matrix(np.array([astral_MTs[-1]]), cell)[0]
#             k = 2
#             result = np.argpartition(dist, k)
#             cell_small1, cell_small2 = result[:k]
#             small_dist1, small_dist2 = dist[result[:k]]
#             # Closest point on a cell
#             astral_MTs[-1] = cell[cell_small1]
#
#             # print("after tang", LA.norm(astral_MTs[-1]-pole))
#         else:
#             astral_MTs[-1], _ = intersect_cell(a, b, get_astral_angle(dummy), pole, cell)
#             # print("after intersect", LA.norm(astral_MTs[-1]-pole))
#         slip = slip - 1 * config.natural_spacing
#         i = i + 1
#         if (i > 20):
#             break
#         # print(f'SLIP FINAL={LA.norm(astral_MTs[-1]-pole)}, coords={astral_MTs[-1]}')
#     return astral_MTs


def smoothen_cell(cell):
    # Example: Replace this with your actual contour coordinates (N, 2)
    contours = cell  # Replace with actual contour coordinates

    # Separate x and y coordinates
    x = contours[:, 0]
    y = contours[:, 1]

    # Apply Gaussian smoothing to both x and y coordinates
    sigma = 10  # Adjust sigma for more or less smoothing
    x_smooth = gaussian_filter1d(x, sigma=sigma)
    y_smooth = gaussian_filter1d(y, sigma=sigma)

    # Combine the smoothed x and y coordinates back into an array
    smoothed_contours = np.column_stack((x_smooth, y_smooth))

    # Plot the original and smoothed contours for comparison
    # plt.figure(figsize=(6, 6))
    # plt.plot(x, y, label='Original Contours', color='red', linewidth=2)
    # plt.plot(x_smooth, y_smooth, label='Smoothed Contours', color='blue', linewidth=2)
    # plt.legend()
    # plt.gca().set_aspect('equal')
    # plt.title('Original vs Smoothed Contours')
    # plt.show()
    return smoothed_contours


def transfer_astro_which(new_length, b1):
    # Create array a1 of specified shape filled with zeroes

    a1 = np.zeros((new_length, 2))

    # Determine the length of the smaller array
    min_length = min(new_length, len(b1))

    # Transfer elements from b1 to a1
    a1[:min_length] = b1[:min_length]

    return a1


def transfer_free_spots(new_length, b1):
    # Create array a1 of specified shape filled with zeroes

    a1 = np.zeros(new_length)

    # Determine the length of the smaller array
    min_length = min(new_length, len(b1))

    # Transfer elements from b1 to a1
    a1[:min_length] = b1[:min_length]

    return a1


def transform_junctions():
    # Get spindle size from its contours
    if (c > 12):
        spindle_contour = get_contours(
            os.path.join(DATA_DIR, 'spindles', 'spindle_' + str(c), 'Mask_' + str(starts[c - 1]) + '.jpg'))
    else:
        spindle_contour = get_contours(
            os.path.join(DATA_DIR, 'spindles', 'spindle_' + str(c), 'Spindle_' + str(starts[c - 1]) + '.jpg'))
    spindle_input = np.zeros((len(spindle_contour), 2))
    for i in range(len(spindle_input)):
        spindle_input[i, 0] = spindle_contour[i, 0] / rescale[c - 1]
        spindle_input[i, 1] = spindle_contour[i, 1] / rescale[c - 1]
    shapely_string = LineString(spindle_input)
    spindle_center = np.array([shapely_string.centroid.x, shapely_string.centroid.y])

    # print(f'spindle center -{20.61,31.83} -> {spindle_center}')
    measured_spindle_center = np.array([[20.61, -31.83],
                                        [23.6, -47.2],
                                        [16.15, -52.10],
                                        [0, 0],
                                        [14.59, -38.55],
                                        [0, 0],
                                        [17.87, -25.61],
                                        [26.79, -28.4],
                                        [14.77, -33.83],
                                        [11.2, -24.55],
                                        [22.29, -55.57],
                                        [16.42, -79.26],
                                        [21.76, -42.18],
                                        [21.33, -49.93],
                                        [27.72, -73.78]])
    # FGs
    file_path = os.path.join(DATA_DIR, "junctions", "C" + str(c) + "_Mask_Movie_Junctions.txt")
    data = np.loadtxt(file_path)
    data[:, 1] = -data[:, 1]
    juncs = data.reshape(-1, 2, 2)

    transf_coef = spindle_center / measured_spindle_center[c - 1]

    juncs = juncs * transf_coef - spindle_center

    return juncs


def rotate_points(points, spindle_angle):
    centroid = np.mean(points, axis=0)
    # Create the 2D rotation matrix
    rotation_matrix = np.array([
        [np.cos(spindle_angle), np.sin(spindle_angle)],
        [-np.sin(spindle_angle), np.cos(spindle_angle)]
    ])
    # Translate points to origin (subtract centroid)
    translated_points = points - centroid
    # Apply the rotation matrix
    rotated_points = np.dot(translated_points, rotation_matrix)
    rotated_points += centroid
    return rotated_points


def slice_array(array, start, end, num_points):
    if num_points < 2:
        raise ValueError("num_points must be at least 2 to form a slice with distinct start and end points")

    # Calculate the step size
    step = (end - start) / (num_points - 1)

    # Generate the indices for the slice
    indices = np.linspace(start, end, num_points).astype(int)

    # indices = gaussian_indexes_simple(start, end, num_points)
    # indices=gaussian_spaced_indices(start, end, num_points)
    # Use the indices to slice the array
    return array[indices]


def distance_matrix(points1, points2):
    # points1: (N1, 2), points2: (N2, 2)
    d = points1[:, None, :] - points2[None, :, :]
    return np.linalg.norm(d, axis=2)


def make_astral_MTs(params, cell, spindle_poles, spindle_angle, spots):
    a = params[0][0]
    r = params[0][2]  # spindle length=2*r
    FG_density = int(params[1])
    n_astro = int(params[2])
    b = params[0][1]

    # Initializing all arrays to keep MTs data

    # state=np.ones((2,(n_astro))) #two states: growing=1, shrinking=-1
    np.random.seed(42)
    free_spots = np.zeros((len(spots)))  # array 0 if spot is free, 1 is taken
    astral_which_spot = np.zeros((len(spots),
                                  2)) - 1  # which astro occupies given spot by index (astral_which_spot[5]=[1,20] means the fifth spot is occupied by (1,20)), i dont remember where -1 is coming from, its coming from not confusing value 0,0 (-1,-1) with actual MT [0,0]
    which_bind = np.zeros((2, int(n_astro)))  # astral MTs have binded
    which_push = np.zeros((2, int(n_astro)))  # astral MT that push
    astral_MTs = np.zeros((2, int(n_astro), config.discr, 2))

    # Astral MTs angles
    astral_angles = np.zeros((2, n_astro))
    astral_angles[0] = np.linspace(spindle_angle - spread / 2, spindle_angle + spread / 2, n_astro)  # predefined angles
    astral_angles[1] = np.linspace(spindle_angle + np.pi - spread / 2, spindle_angle + np.pi + spread / 2, n_astro)

    # orig_length=astral_initial_length*np.ones((5,n_astro)) #to keep track of astral MTs lengthes in case of elongation

    # orig_length=np.zeros((2,n_astro))+0.05
    # Generate the random array
    # orig_length = abs(np.random.normal(loc=mean, scale=stdev, size=(2, n_astro)))
    # orig_length=1*np.ones((2,n_astro)) #to keep track of astral MTs lengthes in case of elongation

    # max_length=np.zeros((2,n_astro))
    # for i in range (2):
    #     for j in range (n_astro):
    #         length_tc, _ =intersect_cell(a,b,astral_angles[i,j],spindle_poles[i],cell)
    #         max_length[i,j]=LA.norm(length_tc-spindle_poles[i])
    if (config.length_MTs == 'gamma'):
        orig_length = np.random.gamma(AL, scale, size=(2, n_astro))  # *max_length
    else:
        orig_length = AL * np.ones((2, n_astro))
    # if (config.length_MTs=='gamma'):
    #     orig_length = np.random.gamma(AL, scale, size=(2, n_astro))
    # else:
    #     orig_length=AL*np.ones((2,n_astro))
    total_rate = config.rescue_rate + config.catastr_rate
    p_grow = config.rescue_rate / total_rate
    p_shrink = config.catastr_rate / total_rate
    if (config.state_MTs == 'random'):
        state = np.random.choice([-1, 1], size=(2, n_astro),
                                 p=[p_shrink, p_grow])  # Randomly assign -1 or 1 based on probabilities
    else:
        state = np.ones((2, n_astro))
    df_list2 = []

    for i in range(2):
        for j in range(n_astro):
            # length_tc, _ =intersect_cell(a,b,astral_angles[i,j],spindle_poles[i],cell)
            # orig_length[i,j]=LA.norm(length_tc-spindle_poles[i])*bounded_normal_random(mean, stdev)

            # print(f'astral=[{i,j}]')
            end = grow_astralMT(a, b, astral_angles[i, j], spindle_poles[i], cell, orig_length[i, j])
            # state[i,j]=1 if (np.random.rand()>0.5) else -1
            astral_MTs[i, j, :, 0] = np.linspace(spindle_poles[i, 0], end[0], config.discr)
            astral_MTs[i, j, :, 1] = np.linspace(spindle_poles[i, 1], end[1], config.discr)

            # which_bind[i,j],free_spots,astral_which_spot=check_bind(i, j,astral_MTs[i,j,-1],spindle_poles,spots,free_spots,astral_which_spot) #used to be check_bind_init
            # if (which_bind[i,j]==1):
            #     which_push[i,j]=0
            #     state[i,j]=-1
            # # state[i,j]=1 if which_bind[i,j]==0 else -1
            # which_push[i,j]=check_push(a,b,astral_MTs[i,j,-1],state[i,j],astral_angles[i,j],spindle_poles[i],cell)
            # which_push[i,j]=check_push_junc(a,b,astral_MTs[i,j,-1],which_bind[i,j],state[i,j],astral_angles[i,j],spindle_poles[i],cell, spots)

            which_bind[i, j], which_push[i, j], state[i, j], free_spots, astral_which_spot = check_push_bind(i, j,
                                                                                                             astral_MTs[
                                                                                                                 i, j, -1],
                                                                                                             spindle_poles,
                                                                                                             spots,
                                                                                                             free_spots,
                                                                                                             astral_which_spot,
                                                                                                             state[
                                                                                                                 i, j],
                                                                                                             astral_angles[
                                                                                                                 i, j],
                                                                                                             cell)

            orig_length[i, j] = LA.norm(astral_MTs[i, j, -1] - spindle_poles[i])

            # List to keep MTs data
            df_list = []
            df_list.append(-1)
            df_list.append([i, j])
            df_list.append('Not dead')
            df_list.append('Neither R or C')
            df_list.append('NOT TOO SHORT')
            df_list.append(which_bind[i, j])
            df_list.append(which_push[i, j])
            df_list.append(state[i, j])
            df_list.append(0)
            df_list.append(LA.norm(astral_MTs[i, j, -1] - spindle_poles[i]))
            df_list.append(orig_length[i, j])
            df_list.append(math.sqrt(
                astral_MTs[i, j, -1, 0] * astral_MTs[i, j, -1, 0] + astral_MTs[i, j, -1, 1] * astral_MTs[i, j, -1, 1]))
            df_list.append(math.degrees((astral_angles[i, j])))
            df_list2.append(df_list)

    return astral_MTs, astral_angles, state, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2


# ============================================================
# SIMULATION PARAMETERS
# ============================================================
class SimulationParameters:
    def __init__(self):
        # Model parameters
        self.number_of_sides = 480
        self.length_MTs = 'gamma'  # default astral MT length distribution for all cell types

        self.state_MTs = 'random'

        self.max_step = 0.1
        self.max_interact_dist = 0.02

        self.push_dist = 0.0005  # 1e-15#0.02
        self.min_cortex_dist = 0.01  # 25
        self.MT_min_length = 0.01  # 25
        self.MT_max_length = 10

        self.natural_spacing = 2 * np.pi / self.number_of_sides
        self.elongate_limit = 0.5
        self.repel_dist = 0.005
        self.discr = 2
        self.astro_lin_density = 1 / 5
        self.max_slip = np.pi / 6

        # Experiment parameters
        if (cell_type == 'FE' or cell_type == 'endo'):
            self.catastr_rate = 0.021 * RC_rate
            self.rescue_rate = 0.029 * RC_rate
            self.shrink_rate = -0.027  * GS_rate
            self.growth_rate = 0.013  * GS_rate
        elif (cell_type == 'celegans'):
            self.catastr_rate = 0.014
            self.rescue_rate = 0.044  # 29
            self.shrink_rate = -18 / 600
            self.growth_rate = 9 / 600

        # Motor protein binding/unbinding
        self.dyn_bind = 0.03
        self.dyn_unbind = 0.02

        # Force
        self.pull = pull
        self.push = push
        self.repel = 0
        # Rigidity & Viscosity
        self.EI = 20 * 0.1 * 0.1
        self.visc = 100
        self.mu_fric = 500

        # Features
        self.buckling = 0
        # self.slipping = 0
        self.pivoting = 0
        self.mobile_motors = 0
        self.force_velocity = 0
        self.v_0 = 0.086
        self.show_vectors = 0
        self.uniform_status = 0


    def save_parameters_to_file(self, directory, filename, additional_params=None):
        """
        Save all parameters of the class and additional parameters to a text file.

        Args:
            directory (str): The directory where the file should be saved.
            filename (str): The name of the file to save the parameters to.
            additional_params (dict): Additional parameters to save, as a dictionary.
        """
        # Create the full file path
        filepath = os.path.join(directory, filename)

        # Check if the file already exists
        if os.path.exists(filepath):
            print(f"File '{filepath}' already exists. Skipping save.")
            return

        # Create the directory if it doesn't exist
        os.makedirs(directory, exist_ok=True)

        # Save the parameters to the file
        with open(filepath, 'w') as file:
            # Write class attributes
            for attr, value in self.__dict__.items():
                file.write(f"{attr}: {value}\n")

            # Write additional parameters
            if additional_params:
                for key, value in additional_params.items():
                    file.write(f"{key}: {value}\n")

        print(f"Parameters saved to '{filepath}'.")


# ============================================================
# SMALL SHARED PHYSICS FUNCTIONS WITH CELL-TYPE BRANCH POINTS
# ============================================================
def gauss_points(N):
    # Gaussian parameters from the fit (replace with your actual values)
    if (cell_type == 'celegans'):
        mu = 0
        sigma = 60
    else:
        """Taken from Mud distribution quantification in follicle cells K. Neville et al. 2023 The EMBO reports"""
        mu = 100.7  # Mean of the Gaussian
        sigma = 30.3  # Standard deviation of the Gaussian

    # Generate 100 uniformly spaced points in [0, 1]
    uniform_points = np.linspace(0, 1, N)

    # Map the uniform points to [0, 180] using the inverse CDF (quantile function) of the Gaussian
    gaussian_spaced_points = norm.ppf(uniform_points, loc=mu, scale=sigma)

    # Clip points to ensure they stay within [0, 180] (not applied for celegans, whose FGs span the full cortex)
    if (cell_type != 'celegans'):
        gaussian_spaced_points = np.clip(gaussian_spaced_points, 0, 180)

    # Sort the points (optional, but ensures they are in ascending order)
    gaussian_spaced_points = np.sort(gaussian_spaced_points)
    return gaussian_spaced_points


def chrom(x, y, angle):
    if (cell_type == 'FE'):
        rc = 0.1
    else:
        rc = 0.05

    line1 = np.array([[x - rc * np.cos(angle - np.pi / 4), x + rc * np.cos(angle - np.pi / 4)],
                      [y - rc * np.sin(angle - np.pi / 4), y + rc * np.sin(angle - np.pi / 4)]])
    line2 = np.array([[x - rc * np.cos(angle + np.pi / 4), x + rc * np.cos(angle + np.pi / 4)],
                      [y - rc * np.sin(angle + np.pi / 4), y + rc * np.sin(angle + np.pi / 4)]])
    return line1, line2


def new_spindle_poles(spindle_poles, spindle_angle, V, Omega):
    if (cell_type == 'celegans'):
        # celegans spindle length can change over time (elongation), so re-derive r from the current poles
        r_local = np.linalg.norm(spindle_poles[0] - spindle_poles[1]) / 2  # radius of the spindle
    else:
        r_local = r

    new_poles = spindle_poles + (V * time_step)
    new_angle = spindle_angle + (Omega * time_step)

    com = np.array([(new_poles[0, 0] + new_poles[1, 0]) / 2,
                    (new_poles[0, 1] + new_poles[1, 1]) / 2])  # rotation is about com(centre of mass)

    new_poles[0] = [r_local * np.cos(new_angle) + com[0], r_local * np.sin(new_angle) + com[1]]
    new_poles[1] = [r_local * np.cos(new_angle + np.pi) + com[0], r_local * np.sin(new_angle + np.pi) + com[1]]
    return new_poles, new_angle


def find_torque(spindle_poles, r, force1, force2):
    if (cell_type == 'celegans'):
        # celegans spindle length can change over time (elongation), so re-derive r from the current poles
        r = np.linalg.norm(spindle_poles[0] - spindle_poles[1]) / 2  # radius of the spindle

    poles12_unit = (spindle_poles[0] - spindle_poles[1]) / (2 * r)

    force1_t = force1 - force1.dot(poles12_unit) * poles12_unit  # transverse to the spindle
    force2_t = force2 - force2.dot(-poles12_unit) * (-poles12_unit)

    # NumPy >=2.0 removed support for 2-D vectors in np.cross (it now requires
    # 3-D vectors). These are 2-D (x,y) vectors in the cell plane, so we
    # compute the scalar z-component of the cross product directly, which is
    # exactly what np.cross(a, b) returned for 2-D inputs pre-2.0.
    sign1 = np.sign(
        poles12_unit[0] * force1_t[1] - poles12_unit[1] * force1_t[0])  # sign for torque using sign of cross product
    sign2 = np.sign((-poles12_unit[0]) * force2_t[1] - (-poles12_unit[1]) * force2_t[0])

    arm1, arm2 = r, r

    # Floating point error may cause dot product to go outside [-1,1]

    if (np.linalg.norm(force1) != 0):
        input_value_safe = np.clip(np.dot(poles12_unit, force1 / np.linalg.norm(force1)), -1,
                                   1)  # Clip the value to be within [-1, 1]
        angle1 = np.arccos(input_value_safe)

    if (np.linalg.norm(force2) != 0):
        input_value_safe = np.clip(np.dot(-poles12_unit, force2 / np.linalg.norm(force2)), -1,
                                   1)  # Clip the value to be within [-1, 1]
        angle2 = np.arccos(input_value_safe)

    torque = sign1 * LA.norm(force1_t) * arm1 + sign2 * LA.norm(force2_t) * arm2

    return torque


def find_force(astral_MTs, spindle_poles, which_bind, which_push, v_c):
    # Initialize force vectors
    push_vecs = np.zeros((2, len(astral_MTs[0]), 2))
    pull_vecs = np.zeros((2, len(astral_MTs[0]), 2))

    # Nested loop to go through each MT
    for i in range(2):
        for j in range(len(astral_MTs[0])):

            # Each force vector is collinear with astral MT body alignment
            vec = np.multiply(astral_MTs[i, j, -1] - spindle_poles[i],
                              1 / LA.norm(astral_MTs[i, j, -1] - spindle_poles[i]))

            # Case 1: MT is being pulled on
            if which_bind[i, j] == 1 and which_push[i, j] == 0:

                if config.force_velocity == 1:
                    a_vec = (astral_MTs[i, j, -1] - astral_MTs[i, j, 0]) / LA.norm(
                        astral_MTs[i, j, -1] - astral_MTs[i, j, 0])
                    force = config.pull * (1 - np.dot(v_c, a_vec) / config.v_0)
                else:
                    force = config.pull

                pull_vecs[i, j] = np.multiply(vec, force)

            # Case 2: MT is pushing against the cortex.
            # Zebrafish ('endo') treats any non-bound MT as pushing and uses a length-proportional force;
            # follicle/celegans require which_push==1 and use an Euler-buckling-limited force.
            elif (which_bind[i, j] == 0) if cell_type == 'endo' else (which_bind[i, j] == 0 and which_push[i, j] == 1):

                push = config.push
                astral_len = LA.norm(astral_MTs[i, j, -1] - spindle_poles[i])
                if cell_type == 'endo':
                    force = astral_len * push
                else:
                    force = -min(push, config.EI * np.pi * np.pi / astral_len / astral_len)

                push_vecs[i, j] = np.multiply(vec, force)

    return push_vecs, pull_vecs


def check_push_bind(i, j, astral, spindle_poles, spots, free_spots, astral_which_spot, state, astral_angles, cell):
    """
    Combined function that checks both binding and pushing conditions for an astral MT.
    Reordered workflow: first checks proximity to cell end, then binding distance.

    Parameters:
        i, j: MT indices
        astral: Astral MT tip coordinates
        spindle_poles: Spindle poles coordinates
        spots: FG spot coordinates
        free_spots: Array indicating available spots
        astral_which_spot: Tracking array
        state: MT state (1=growth, -1=shrink)
        astral_angles: MT angles
        cell: Cell geometry

    Returns:
        bind: 1 if bound, 0 otherwise
        push: 1 if pushing, 0 otherwise
        free_spots: Updated free spots array
        astral_which_spot: Updated tracking array
    """
    # Initialize outputs
    bind = 0
    push = 0

    dist = distance_matrix(np.array([astral]), spots)
    # check if the MT is close enough to a FG spot
    # and if the spot is free
    if (np.min(dist[0]) <= config.max_interact_dist and free_spots[np.argmin(dist[0])] == 0 and random.uniform(0,
                                                                                                               1) <= prob_dyn_bind):

        free_spots[np.argmin(dist[0])] = 1
        bind = 1
        state = -1
        astral_which_spot[np.argmin(dist[0]), 0] = i
        astral_which_spot[np.argmin(dist[0]), 1] = j
    else:
        end, _ = intersect_cell(a, b, astral_angles, spindle_poles[i], cell)
        # Check if the astral MT is close enough to the cortex.
        # celegans uses the magnitude of the difference vector; follicle/zebrafish use the difference of magnitudes.
        if (cell_type == 'celegans'):
            close_enough = math.isclose(abs(LA.norm(end - astral)), 0, abs_tol=config.push_dist)
        else:
            close_enough = math.isclose(abs(LA.norm(end) - LA.norm(astral)), 0, abs_tol=config.push_dist)
        if close_enough:
            # If not binding, check if the MT is pushing
            push = 1
            state = 1

    return bind, push, state, free_spots, astral_which_spot


def slide_and_rotate_spindle(spindle_poles, spindle_angle, V, Omega, cell, r, w, config, dt):
    """
    Attempts maximal safe translation and then maximal safe rotation of the spindle in the cell.
    Translation and rotation are decoupled: if translation is impossible, still attempt full rotation.
    Returns new poles, angle, and a boolean (True if any move occurred, else False).
    """
    print("\n[DEBUG] Starting slide_and_rotate_spindle")
    print(f"[DEBUG] Input spindle_angle: {np.rad2deg(spindle_angle):.2f} deg, V: {V}, Omega: {Omega}, dt: {dt}")

    # Zebrafish uses fewer bisection iterations (10) than follicle/celegans (20)
    n_iter = 10 if cell_type == 'endo' else 20

    # --- Maximal translation in the direction of V (not just tangential) ---
    v_norm = np.linalg.norm(V)
    step_size = min(v_norm * dt, config.max_step) if v_norm > 0 else 0
    print(f"[DEBUG] Velocity: {V}, norm: {v_norm}, step_size: {step_size}")

    low, high = 0, 1
    for j in range(n_iter):
        frac = (low + high) / 2
        trial_poles = spindle_poles + frac * V / v_norm * step_size if v_norm > 0 else spindle_poles.copy()
        is_valid = check_spindle(trial_poles, spindle_angle, cell, r, w)
        if is_valid:
            low = frac
        else:
            high = frac
        print(f"[DEBUG] Trans iter {j}: frac={frac:.4f}, valid={is_valid}")
    spindle_poles_trans = spindle_poles + low * V / v_norm * step_size if v_norm > 0 else spindle_poles.copy()
    moved_translation = np.any(np.abs(spindle_poles_trans - spindle_poles) > 1e-9)

    # --- Maximal allowed rotation at the translated position ---
    dtheta = Omega * dt
    low_rot, high_rot = 0, 1
    for k in range(n_iter):
        frac = (low_rot + high_rot) / 2
        trial_angle = spindle_angle + frac * dtheta
        com = np.mean(spindle_poles_trans, axis=0)
        R = np.array([
            [np.cos(trial_angle - spindle_angle), -np.sin(trial_angle - spindle_angle)],
            [np.sin(trial_angle - spindle_angle), np.cos(trial_angle - spindle_angle)]
        ])
        trial_poles_rot = (spindle_poles_trans - com) @ R.T + com
        is_valid = check_spindle(trial_poles_rot, trial_angle, cell, r, w)
        if is_valid:
            low_rot = frac
        else:
            high_rot = frac
        print(f"[DEBUG] Rot iter {k}: frac={frac:.4f}, angle={np.rad2deg(trial_angle):.2f} deg, valid={is_valid}")

    final_angle = spindle_angle + low_rot * dtheta
    com = np.mean(spindle_poles_trans, axis=0)
    R = np.array([
        [np.cos(final_angle - spindle_angle), -np.sin(final_angle - spindle_angle)],
        [np.sin(final_angle - spindle_angle), np.cos(final_angle - spindle_angle)]
    ])
    final_poles = (spindle_poles_trans - com) @ R.T + com
    moved_rotation = (abs(final_angle - spindle_angle) > 1e-9)

    # --- Determine if any change occurred and return ---
    has_changed = (moved_translation or moved_rotation)

    if not check_spindle(final_poles, final_angle, cell, r, w):
        print("[WARNING] No valid translation/rotation possible; returning original spindle configuration.")
        return spindle_poles.copy(), spindle_angle, False
    elif has_changed:
        print("[DEBUG] --- End slide_and_rotate_spindle: MOVE occurred ---\n")
        return final_poles, final_angle, True
    else:
        print("[DEBUG] --- End slide_and_rotate_spindle: NONE ---\n")
        return spindle_poles.copy(), spindle_angle, False


# ============================================================
# ASTRAL MICROTUBULE UPDATE (shared, with cell-type branch points)
# ============================================================
def update_astral_MTs(params, cell, spindle_poles, spindle_angle, delta_spindle_angle, astral_MTs, astral_angles, state,
                      which_push, which_bind, spots, free_spots, astral_which_spot, orig_length, force_vector_1,
                      force_vector_2, run):
    a = params[0][0]
    r = params[0][2]  # spindle length=2*r
    FG_density = int(params[1])
    n_astro = int(params[2])
    b = params[0][1]
    t_time = run * time_step

    # Updating astral MTs angles
    """
    Pivoting is the process of changing the angle of the astral MTs from its original angle.
    If pivoting is on, astral microtubule will stay at the same angle.
    If pivoting is off, the astral MTs will return to their original angle.
    """

    """
    Refreshing astral MTs angles for a new time step. Used as a placeholder for the new angles.
    """
    if (config.pivoting == 0):
        # Astral MTs angles are updated to the new spindle angle and they are constant in reference to the spindle.
        astral_angles[0] = np.linspace(spindle_angle - spread / 2, spindle_angle + spread / 2, n_astro)
        astral_angles[1] = np.linspace(spindle_angle + np.pi - spread / 2, spindle_angle + np.pi + spread / 2, n_astro)
    else:
        # Astral MTs angles change together with the spindle angle.
        astral_angles[0] = astral_angles[0] + delta_spindle_angle
        astral_angles[1] = astral_angles[1] + delta_spindle_angle

    # List to keep MTs data
    df_list2 = []

    for i in range(2):
        for j in range(n_astro):
            df_list = []  # append 'Number','Death', 'Switch','short','Push','State','Length','end if out','astral angle'
            df_list.append(run)
            df_list.append([i, j])
            astral_MTs[i, j, 0] = spindle_poles[
                i]  # Updating MT start, but this stretches microtubule because the plus end is not updated
            # Ways MT length can change:
            # 1. pulling MTs stretched
            # 2. growing MTs +end grow
            # 3. shrinking MTs -end shrink

            # Check if any has unbinded

            if (which_bind[i, j] == 1):

                astral_MTs[i, j, 0] = spindle_poles[
                    i]  # Updating MT start, but this stretches microtubule because the plus end is not updated

                if (config.pivoting == 1):
                    """Astral MTs angles are updated to the new spindle angle after spindle has moved (and microtubule minus end correspondingly)."""
                    astral_angles[i, j] = get_astral_angle(astral_MTs[i, j])  # config.pivoting allowed
                astral_intersect, geom_type = intersect_cell(a, b, astral_angles[i, j], spindle_poles[i], cell)

                # Check if
                if (random.uniform(0, 1) <= prob_dyn_unbind or geom_type == 'MultiPoint'):
                    which_bind[i, j] = 0
                    free_spots[np.argwhere((astral_which_spot[:, 0] == i) & (
                                astral_which_spot[:, 1] == j))] = 0  # book keeping for available FG slots
                    state[i, j] = -1
                    if (geom_type == 'MultiPoint'):
                        astral_MTs[i, j, -1] = astral_intersect
                        astral_MTs[i, j] = restructure(astral_MTs[i, j])
                    df_list.append('Unbinded')
                    df_list.append('Neither R or C')
                    df_list.append('unbinded')
                else:
                    if (cell_type == 'endo'):
                        # Zebrafish: as the tracked cell shape/FG layout is refreshed periodically, keep the
                        # bound astral MT glued to (or re-find) its FG spot.
                        if (
                                config.mobile_motors == 1):  # if the cell shape changes and new spots forms keep the astral MT glued to the same index spot
                            astral_MTs[i, j, -1] = spots[
                                np.where((astral_which_spot == np.array([i, j])).all(axis=1))[0]]

                        else:
                            free_spots[np.argwhere((astral_which_spot[:, 0] == i) & (
                                        astral_which_spot[:, 1] == j))] = 0  # release old spot

                            dist = distance_matrix(np.array([astral_intersect]), spots)  # find new closest spot

                            astral_MTs[i, j, -1] = spots[np.argmin(dist[0])]  # put the end to the new spot
                            which_bind[i, j], which_push[i, j], state[
                                i, j], free_spots, astral_which_spot = check_push_bind_init(i, j, astral_MTs[i, j, -1],
                                                                                            spindle_poles, spots,
                                                                                            free_spots,
                                                                                            astral_which_spot,
                                                                                            state[i, j],
                                                                                            astral_angles[i, j], cell)

                        astral_MTs[i, j] = restructure(astral_MTs[i, j])

                    df_list.append('Not Dead not binded')
                    df_list.append('Neither R or C')
                    df_list.append('stays binded')

            # NOT BINDED
            elif (which_bind[i, j] == 0):  # didn't bind (which_bind[i,j]=0)
                df_list.append('Not Dead not binded')
                rand_n = random.uniform(0, 1)
                if (state[i, j] > 0 and rand_n <= prob_catastr):  # if catastrophe happens
                    state[i, j] = -1
                    rate = config.shrink_rate
                    df_list.append('Catastrophe')
                elif (state[i, j] < 0 and rand_n <= prob_rescue):  # if rescue happens
                    state[i, j] = 1
                    rate = config.growth_rate
                    df_list.append('Rescue')
                else:
                    df_list.append('Neither R or C')
                    rate = config.shrink_rate if state[i, j] == -1 else config.growth_rate

                # GROWING AND SHRINKING

                # Floating shrinking and growing MTs spin with the spindle

                # SHRINKING
                if (state[i, j] == -1):
                    old_astral_MTs = spindle_poles[i] + orig_length[i, j] * np.array([np.cos(astral_angles[i, j]),
                                                                                      np.sin(astral_angles[
                                                                                                 i, j])])  # Shifting astral MTs due to new -end coordinate
                    length_delta = time_step * rate * np.array(
                        [np.cos(astral_angles[i, j]), np.sin(astral_angles[i, j])])  # shrink length
                    astral_MTs[i, j, -1] = old_astral_MTs + length_delta  # Substracting
                    astral_MTs[i, j, 0] = spindle_poles[
                        i]  # Updating MT start, but this stretches microtubule because the plus end is not updated
                    astral_MTs[i, j] = restructure(astral_MTs[i, j])
                    astral_intersect, _ = intersect_cell(a, b, astral_angles[i, j], spindle_poles[i], cell)
                    # Comparing which point is further from the spindle pole (intersect or shrinking end)
                    if (LA.norm(astral_MTs[i, j, -1] - spindle_poles[i]) > LA.norm(
                            astral_intersect - spindle_poles[i])):
                        astral_MTs[i, j, -1] = astral_intersect
                        astral_MTs[i, j] = restructure(astral_MTs[i, j])

                    if (LA.norm(astral_MTs[i, j, -1] - spindle_poles[
                        i]) < config.MT_min_length):  # Check if astral MT is too short and is going to be replaced
                        # Angle Options for new astral
                        opt_astral_angles = np.zeros((2, n_astro))
                        opt_astral_angles[0] = np.linspace(spindle_angle - spread / 2, spindle_angle + spread / 2,
                                                           n_astro)
                        opt_astral_angles[1] = np.linspace(spindle_angle + np.pi - spread / 2,
                                                           spindle_angle + np.pi + spread / 2, n_astro)

                        astral_MTs[i, j] = make_new_MT(i, j, a, b, astral_MTs[i, j], spindle_poles, opt_astral_angles,
                                                       cell, free_spots, astral_which_spot, orig_length[i, j])
                        astral_angles[i, j] = get_astral_angle(astral_MTs[i, j])
                        orig_length[i, j] == LA.norm(astral_MTs[i, j, -1] - spindle_poles[i])
                        state[i, j] = 1
                        df_list.append('short and reborn')
                    else:
                        df_list.append('not short shrinking')
                    # Update length and check for a new pulling or pushing MT emergence
                    orig_length[i, j] = LA.norm(astral_MTs[i, j, -1] - spindle_poles[i])

                    which_bind[i, j], which_push[i, j], state[i, j], free_spots, astral_which_spot = check_push_bind(i,
                                                                                                                     j,
                                                                                                                     astral_MTs[
                                                                                                                         i, j, -1],
                                                                                                                     spindle_poles,
                                                                                                                     spots,
                                                                                                                     free_spots,
                                                                                                                     astral_which_spot,
                                                                                                                     state[
                                                                                                                         i, j],
                                                                                                                     astral_angles[
                                                                                                                         i, j],
                                                                                                                     cell)

                # GROWING
                else:
                    # old_angle = astral_angles[i, j]
                    # virtual_astral = grow_astralMT(a, b, astral_angles[i, j], spindle_poles[i], cell, orig_length[
                    #     i, j] + config.growth_rate * time_step)  # growing astral MTs without considering config.slipping behaviour
                    # _, virtual_push, _, _, _ = check_push_bind(i, j, astral_MTs[i, j, -1], spindle_poles, spots,
                    #                                            free_spots, astral_which_spot, state[i, j],
                    #                                            astral_angles[i, j], cell)

                    # if (
                    #         config.slipping == 1 and virtual_push == 1):  # Figure out if the MT is touching the cortex before growing (prerequisite for config.slipping)
                    #
                    #     virtual_astral = grow_astralMT(a, b, astral_angles[i, j], spindle_poles[i], cell,
                    #                                    orig_length[i, j])
                    #     center = np.array([1 / 2 * (spindle_poles[0, 0] + spindle_poles[1, 0]),
                    #                        1 / 2 * (spindle_poles[0, 1] + spindle_poles[1, 1])])
                    #     spin_vec = center - spindle_poles[i]
                    #     a_vec = astral_MTs[i, j, -1] - astral_MTs[i, j, 0]
                    #     as_angle = angle_between_vectors(a_vec, spin_vec)
                    #
                    #     cut_part = max(0, (orig_length[i, j] + config.growth_rate * time_step) - LA.norm(
                    #         virtual_astral - spindle_poles[i]))
                    #     astral_MTs[i, j] = slipping_f(i, j, a, b, cell, astral_MTs[i, j], spindle_poles[i],
                    #                                   orig_length[i, j], spindle_angle, n_astro,
                    #                                   cut_part)  # cell, astral_MTs,pole,orig_length
                    #     astral_angles[i, j] = get_astral_angle(astral_MTs[i, j])
                    #     which_bind[i, j], which_push[i, j], state[
                    #         i, j], free_spots, astral_which_spot = check_push_bind(i, j, astral_MTs[i, j, -1],
                    #                                                                spindle_poles, spots, free_spots,
                    #                                                                astral_which_spot, state[i, j],
                    #                                                                astral_angles[i, j], cell)
                    #
                    #     df_list.append('not short config.slipping')
                    # else:
                    astral_MTs[i, j, -1] = grow_astralMT(a, b, astral_angles[i, j], spindle_poles[i], cell,
                                                         orig_length[i, j])
                    which_bind[i, j], which_push[i, j], state[
                        i, j], free_spots, astral_which_spot = check_push_bind(i, j, astral_MTs[i, j, -1],
                                                                               spindle_poles, spots, free_spots,
                                                                               astral_which_spot, state[i, j],
                                                                               astral_angles[i, j], cell)
                    df_list.append('not short growing')
            orig_length[i, j] = LA.norm(astral_MTs[i, j, -1] - spindle_poles[i])

            astral_intersect, _ = intersect_cell(a, b, astral_angles[i, j], spindle_poles[i], cell)

            # Recording data to excel sheet
            df_list.append(which_bind[i, j])
            df_list.append(which_push[i, j])
            df_list.append(state[i, j])
            if (i == 0):
                df_list.append(np.linalg.norm(force_vector_1[j]))
            else:
                df_list.append(np.linalg.norm(force_vector_2[j]))

            df_list.append(
                LA.norm([astral_MTs[i, j, -1, 0] - spindle_poles[i, 0], astral_MTs[i, j, -1, 1] - spindle_poles[i, 1]]))
            df_list.append(orig_length[i, j])
            df_list.append(math.sqrt(
                astral_MTs[i, j, -1, 0] * astral_MTs[i, j, -1, 0] + astral_MTs[i, j, -1, 1] * astral_MTs[i, j, -1, 1]))
            df_list.append(math.degrees((astral_angles[i, j])))
            df_list2.append(df_list)
    return astral_MTs, astral_angles, state, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2


# ============================================================
# SPINDLE MOVEMENT (shared, with cell-type branch points)
# ============================================================
def move_spindle(params, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, cell, spots, which_push,
                 which_bind, free_spots, astral_which_spot, orig_length, v_c, i):
    # Parameters
    a = params[0][0]
    b = params[0][1]
    r = params[0][2]  # spindle length=2*r
    n_astro = int(params[2])
    t_time = i * time_step
    og_spindle_angle = spindle_angle

    dum_spot = astral_which_spot
    if (cell_type == 'endo'):
        if (i % (1 / time_step) == 0 and i < total_time / time_step):
            # Zebrafish cell shape is tracked from live movies and periodically refreshed.
            cell = update_cell_zebrafish(t_time)
            if model == 'ever':
                spots = make_fgs_zebrafish_uf(int(params[1]), a, b, cell, spindle_poles,
                                              int(t_time / frame_rates[c - 1]))
            elif model == 'junc':
                spots = make_fgs_zebrafish_junc(int(params[1]), a, b, cell, spindle_poles,
                                                int(t_time / frame_rates[c - 1]))
            free_spots = transfer_free_spots(len(spots), free_spots)
            astral_which_spot = transfer_astro_which(len(spots), dum_spot)

    if (cell_type == 'celegans'):
        spindle_length = LA.norm(spindle_poles[0] - spindle_poles[1])  # Spindle length

    if (cell_type == 'endo'):
        if not check_spindle(spindle_poles, spindle_angle, cell, r, w):
            print("[WARNING] Spindle is outside the cell after adjustment.")
            spindle_poles, spindle_angle = project_spindle_to_feasible(spindle_poles, spindle_angle, cell, r, w)
            delta_spindle_angle = spindle_angle - og_spindle_angle

    push_force, pull_force = find_force(astral_MTs, spindle_poles, which_bind, which_push, v_c)

    # Find TOTAL pulling and pushing forces
    pull_t = np.sum(pull_force[0] + pull_force[1], axis=0)
    push_t = np.sum(push_force[0] + push_force[1], axis=0)

    # Calculating ratio of pulling and pushing forces

    if ((LA.norm(pull_t) + LA.norm(push_t)) == 0):
        ratio = 0
    else:
        ratio = 100 * LA.norm(pull_t) / (LA.norm(pull_t) + LA.norm(push_t))

    force_vector_1 = pull_force[0] + push_force[0]
    force_vector_2 = pull_force[1] + push_force[1]

    # Adding repel force
    force1 = np.sum(force_vector_1, axis=0)
    force2 = np.sum(force_vector_2, axis=0)

    if (cell_type == 'celegans' and t_time < stall_time):
        force_net = 0 * (force1 + force2)
    else:
        force_net = (force1 + force2)

    spindle_dir = np.array([np.cos(spindle_angle), np.sin(spindle_angle)])
    delta_spindle_angle = 0
    mu = config.visc
    torque = find_torque(spindle_poles, r, force1, force2)

    if (cell_type == 'celegans' and model != 'spindle'):
        V = force_net / (6 * np.pi * mu * r)  # Translational velocity for circular cross-section
        Omega = torque / (8 * np.pi * mu * (1 * r) ** 3)
    else:
        V = force_net / (6 * np.pi * mu * r)  # Translational velocity for circular cross-section
        Omega = torque / (8 * np.pi * mu * r ** 3)  #

    old_poles = spindle_poles.copy()
    old_angle = spindle_angle
    if (np.linalg.norm(force_net) != 0 and torque != 0):

        new_poles, new_angle = new_spindle_poles(spindle_poles, spindle_angle, V, Omega)

        if not check_spindle(new_poles, new_angle, cell, r, w):
            print(f"[DEBUG] Spindle INVALID---------------------------------.")
            spindle_poles, spindle_angle, success = slide_and_rotate_spindle(new_poles, new_angle, V, Omega, cell, r, w,
                                                                             config, time_step)
            if not success:
                print(f"[DEBUG] Spindle still INVALID after sliding and rotating. Keeping original poles and angle.")
                spindle_poles, spindle_angle = old_poles, old_angle
        else:
            print(f"[DEBUG] Spindle VALID+++++++++++++++++++++++++++++++++.")
            spindle_poles, spindle_angle = new_poles, new_angle
        delta_spindle_angle = new_angle - spindle_angle

    print(f'[DEBUG] Angle change={np.rad2deg(delta_spindle_angle)}')
    print(f'[DEBUG] Position change={spindle_poles - old_poles, LA.norm(spindle_poles - old_poles)}')

    v_c = (spindle_poles - old_poles) / time_step

    astral_MTs, astral_angles, state, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2 = update_astral_MTs(
        params, cell, spindle_poles, spindle_angle, delta_spindle_angle, astral_MTs, astral_angles, state, which_push,
        which_bind, spots, free_spots, astral_which_spot, orig_length, force_vector_1, force_vector_2, i)

    return cell, spots, spindle_poles, spindle_angle, astral_MTs, astral_angles, state, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2, ratio, push_force, pull_force, v_c


def move_severed_spindle(params, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, cell, spots,
                         which_push, which_bind, free_spots, astral_which_spot, orig_length, v_c, i):
    # Parameters
    a = params[0][0]  # Sphere radius (formerly semi-major axis)
    b = params[0][1]  # Unused now (was semi-minor axis)
    r = params[0][2]  # Original spindle length (now used for distance checks)
    t_time = i * time_step

    # Each centrosome is now an independent sphere with radius 'a'
    mu = config.visc

    # Calculate forces for each centrosome independently
    push_force, pull_force = find_force(astral_MTs, spindle_poles, which_bind, which_push, v_c)

    # Force vectors for each centrosome
    force1 = np.sum(pull_force[0] + push_force[0], axis=0)
    force2 = np.sum(pull_force[1] + push_force[1], axis=0)

    # Calculate movement for each centrosome independently (Stokes' law for spheres)
    def move_sphere(position, force, cell_boundary, sphere_radius):
        # Simple Stokes' drag for sphere
        velocity = force / (6 * np.pi * mu * sphere_radius)

        # Proposed new position
        new_position = position + velocity * time_step

        # Check boundary conditions (similar to original check_spindle but for single sphere)
        if LA.norm(new_position) < cell_boundary - sphere_radius:
            return new_position, velocity
        else:
            # If would hit boundary, scale down movement
            max_dist = cell_boundary - sphere_radius - LA.norm(position)
            if max_dist > 0:
                scale = max_dist / LA.norm(velocity * time_step)
                return position + velocity * time_step * scale, velocity * scale
            else:
                return position, np.zeros(2)

    # Move each centrosome independently. Celegans (post-severance) uses a larger step scale.
    step = 0.25 if cell_type == 'celegans' else 0.025
    cell_radius = LA.norm(cell[0])  # Assuming cell is circular
    new_pole1, v1 = move_sphere(spindle_poles[0], force1, cell_radius, step)
    new_pole2, v2 = move_sphere(spindle_poles[1], force2, cell_radius, step)

    # Update spindle poles (now independent)
    spindle_poles = np.array([new_pole1, new_pole2])
    v_c = np.array([0, 0])  # (spindle_poles - old_poles) / time_step

    # Update spindle angle based on new pole positions
    new_angle = spindle_angle
    delta_spindle_angle = new_angle - spindle_angle
    spindle_angle = new_angle

    # Update MTs (same as before but now poles may be completely independent)
    astral_MTs, astral_angles, state, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2 = update_astral_MTs(
        params, cell, spindle_poles, spindle_angle, delta_spindle_angle, astral_MTs, astral_angles,
        state, which_push, which_bind, spots, free_spots, astral_which_spot, orig_length,
        pull_force[0] + push_force[0], pull_force[1] + push_force[1], i)

    # Calculate force ratios (for reporting)
    pull_t = np.sum(pull_force[0]) + np.sum(pull_force[1])
    push_t = np.sum(push_force[0]) + np.sum(push_force[1])
    ratio = 100 * LA.norm(pull_t) / (LA.norm(pull_t) + LA.norm(push_t)) if (LA.norm(pull_t) + LA.norm(
        push_t)) > 0 else 0

    return cell, spots, spindle_poles, spindle_angle, astral_MTs, astral_angles, state, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2, ratio, push_force, pull_force, v_c


# ============================================================
# Zebrafish-specific helpers: cell shape is tracked from live movies
# (real microscopy contours) rather than generated parametrically.
# ============================================================

def resample_cell(traj, n_points=2400):
    distances = np.cumsum(np.sqrt(np.sum(np.diff(traj, axis=0) ** 2, axis=1)))
    distances = np.insert(distances, 0, 0)  # Add starting point

    # Create interpolation functions for x and y
    fx = interp1d(distances, traj[:, 0], kind='linear')
    fy = interp1d(distances, traj[:, 1], kind='linear')

    # Generate new equally spaced distances
    new_distances = np.linspace(0, distances[-1], n_points)

    # Interpolate new points
    new_x = fx(new_distances)
    new_y = fy(new_distances)

    contour = np.column_stack((new_x, new_y))

    target_point, _ = intersect_cell(a, b, 0, np.array([0, 0]), contour)
    distances = np.linalg.norm(contour - target_point, axis=1)

    # Find the closest point's index
    closest_idx = np.argmin(distances)

    # Roll the array to start at closest_idx
    shifted_contour = np.roll(contour, -closest_idx, axis=0)
    return shifted_contour


def get_real_cell(t_time=0):
    """
    Updating cell shape according to live movies
    """
    # Read the image
    image_path = os.path.join(DATA_DIR, 'cells', 'cell_' + str(c),
                              'Mask_' + str(int(starts[c - 1] + t_time / frame_rates[c - 1])) + '.jpg')
    # Get the raw contours and rescale them
    all_points = get_contours(image_path)

    cell_input = np.zeros((len(all_points), 2))
    for i in range(len(cell_input)):
        cell_input[i, 0] = all_points[i, 0] / rescale[c - 1]
        cell_input[i, 1] = all_points[i, 1] / rescale[c - 1]
    cell = enhance_cell(cell_input)

    # Get spindle size from its contours
    if (c > 12):
        spindle_contour = get_contours(
            os.path.join(DATA_DIR, 'spindles', 'spindle_' + str(c), 'Mask_' + str(starts[c - 1]) + '.jpg'))
    else:
        spindle_contour = get_contours(
            os.path.join(DATA_DIR, 'spindles', 'spindle_' + str(c), 'Spindle_' + str(starts[c - 1]) + '.jpg'))
    spindle_input = np.zeros((len(spindle_contour), 2))
    for i in range(len(spindle_input)):
        spindle_input[i, 0] = spindle_contour[i, 0] / rescale[c - 1]
        spindle_input[i, 1] = spindle_contour[i, 1] / rescale[c - 1]
    spindle_length = max(spindle_input[:, 1]) - min(spindle_input[:, 1])
    r = spindle_length / 2
    shapely_string = LineString(spindle_input)
    spindle_center = np.array([shapely_string.centroid.x, shapely_string.centroid.y])

    cell = cell - spindle_center
    cell = smoothen_cell(cell)
    cell = close_contours(cell)

    cell = resample_cell(cell, n_points=2400)

    return cell


def interpolate_contours(contour1, contour2, steps=5):
    """Generate interpolated contours between two shapes with equal points.

    Args:
        contour1: First contour (N,2 array)
        contour2: Second contour (N,2 array)
        steps: Number of interpolation steps (including endpoints)

    Returns:
        List of interpolated contours
    """
    # Linear interpolation between corresponding points
    alphas = np.linspace(0, 1, steps)
    return [contour1 * (1 - alpha) + contour2 * alpha for alpha in alphas]


def update_cell_zebrafish(t_time):
    steps = frame_rates[c - 1]
    print(f'cell before')
    cell_before = get_real_cell(t_time)
    print(f'cell after')
    cell_after = get_real_cell(t_time + frame_rates[c - 1])

    interpolated = interpolate_contours(cell_before, cell_after, steps)

    current = int((t_time / frame_rates[c - 1] % 1) * steps)
    print(f'current {current}')
    return interpolated[current]


def project_spindle_to_feasible(
        spindle_poles, spindle_angle, cell, r, w,
        step=0.02,  # radial step size (um) for search rings
        max_trans_only=0.1,  # <= THIS caps Stage 1 translation magnitude
        rings=20,  # number of rings to try in Stage 2
        directions=16,  # angular samples per ring
        rot_steps=10,  # rotation samples on each side
        rot_span_deg=20  # total ± span for rotation search (deg)
):
    """
    Project spindle to nearest feasible configuration inside the cell.

    Stage 1: translation-only search, limited to <= max_trans_only displacement.
    Stage 2: translation + rotation search if Stage 1 fails.

    Assumes:
      - generate_spindle(spindle_poles, spindle_angle, r, w) -> (N,2) points
      - config.min_cortex_dist is defined
      - 'cell' is an (M,2) polygon (counterclockwise recommended)
    """
    min_dist = config.min_cortex_dist
    cell_poly = Polygon(cell)

    def is_feasible(poles, angle):
        spindle_body = generate_spindle(poles, angle, r, w)
        spindle_poly = Polygon(spindle_body)
        inside = cell_poly.buffer(-min_dist).contains(spindle_poly)
        sep_ok = (cell_poly.boundary.distance(spindle_poly) >= min_dist)
        return inside and sep_ok

    # ---------------- Stage 1: translation-only (<= max_trans_only) ----------------
    print(f"[DEBUG] Stage 1: Translation-only (limit {max_trans_only} um), step={step}, directions={directions}")
    max_rings_trans_only = int(max_trans_only / step + 1e-9)
    for ring in range(max_rings_trans_only + 1):
        radius = ring * step
        for d in range(directions):
            theta = 2 * np.pi * d / directions
            dx = radius * np.cos(theta)
            dy = radius * np.sin(theta)
            candidate_poles = spindle_poles + np.array([dx, dy])
            feasible = is_feasible(candidate_poles, spindle_angle)
            print(f"[DEBUG][T-only] ring={ring} rad={radius:.3f} dx={dx:.3f} dy={dy:.3f} -> feasible={feasible}")
            if feasible:
                print(f"[DEBUG] Found feasible translation-only move at ring={ring}, dx={dx:.3f}, dy={dy:.3f}")
                return candidate_poles, spindle_angle

    # ---------------- Stage 2: translation + rotation ----------------
    print(f"[DEBUG] Stage 2: Translation + rotation (rings={rings}, rot_steps={rot_steps}, span=±{rot_span_deg}°)")
    for ring in range(rings + 1):
        radius = ring * step
        for d in range(directions):
            theta = 2 * np.pi * d / directions
            dx = radius * np.cos(theta)
            dy = radius * np.sin(theta)
            for rs in range(-rot_steps, rot_steps + 1):
                delta_angle = np.deg2rad(rot_span_deg) * (rs / rot_steps)
                candidate_poles = spindle_poles + np.array([dx, dy])
                candidate_angle = spindle_angle + delta_angle
                feasible = is_feasible(candidate_poles, candidate_angle)
                print(
                    f"[DEBUG][T+R] ring={ring} rad={radius:.3f} dx={dx:.3f} dy={dy:.3f} dθ={np.rad2deg(delta_angle):.2f} -> feasible={feasible}")
                if feasible:
                    print(
                        f"[DEBUG] Found feasible translation+rotation at ring={ring}, dx={dx:.3f}, dy={dy:.3f}, Δθ={np.rad2deg(delta_angle):.2f}°")
                    return candidate_poles, candidate_angle

    print("[DEBUG] No feasible configuration found; returning original.")
    return spindle_poles, spindle_angle


# ============================================================
# FG (motor/spot) placement — cell-type / mode specific
# ============================================================

def make_fgs_follicle(motors, a, b, cell, spindle_poles, frame):
    """
    Generate FG spots for the follicular-epithelial cell (basolateral band,
    mirrored left/right; no apical FGs).
    """
    n_points = motors

    spots = np.zeros((n_points // 2,
                      2))  # make right half of the points first, then reflect and concatenate both sides into one array

    N = n_points // 2
    gaussian_spaced_points = gauss_points(N)
    angles = np.deg2rad(gaussian_spaced_points) - np.pi / 2
    low_limit = np.deg2rad(-60)
    high_limit = np.deg2rad(60)
    angles = angles[(angles >= low_limit) & (angles <= high_limit)]
    # Initialize spots array
    spots = np.zeros((len(angles), 2))

    # Compute x and y coordinates
    spots[:, 0] = a * np.cos(angles)  # x-coordinates
    spots[:, 1] = b * np.sin(angles)  # y-coordinates

    # For squeezed cells to readjust FGs positions because they are based on angles need to readjust angles
    if (a != b):
        for i in range(len(angles)):
            dumb = delta_theta(a, b, angles[i])
            angles[i] = dumb

    # Calculate positions of the spots on the right half of the cell and mirror them to create the left half and concatenate two halves
    spots[:, 0] = a * np.cos(angles - 0 * np.pi)
    spots[:, 1] = b * np.sin(angles - 0 * np.pi)

    other_half = np.zeros((len(angles), 2))  # spots on left side
    other_half[:, 0] = -spots[:, 0]
    other_half[:, 1] = spots[:, 1]
    spots = np.concatenate((spots[:-1], other_half))

    return spots


def make_fgs_neuroblast(motors, a, b, cell, spindle_poles, frame):
    """
    Generate FG spots for the neuroblast (follicular-epithelial geometry with
    NB flag set): apical/top-only FGs, replacing the basolateral band.
    """
    n_points = motors

    # Apical/top-only FGs (the basolateral band computed here in the original
    # script before this NB branch is discarded unused, so it is omitted)
    n_top = n_points
    gaussian_spaced_points = gauss_points(n_top)
    angles = np.deg2rad(gaussian_spaced_points)
    low_limit = np.deg2rad(30)
    high_limit = np.deg2rad(150)
    angles = angles[(angles >= low_limit) & (angles <= high_limit)]

    spots_apical = np.zeros((len(angles), 2))
    spots_apical[:, 0] = a * np.cos(angles)
    spots_apical[:, 1] = b * np.sin(angles)
    return spots_apical


def make_fgs_celegans_pnc(motors, a, b, cell, spindle_poles, frame):
    """
    Generate FG spots for C. elegans pronuclear centering (PNC): FGs spread
    broadly (85-95deg cutoffs) around anterior and posterior poles of the superellipse.
    """
    n_points = motors

    # Anterior FGs
    gaussian_spaced_points = gauss_points(int(0.6 * n_points))
    angles = np.deg2rad(gaussian_spaced_points) + np.pi
    low_limit = np.deg2rad(85)
    high_limit = 2 * np.pi - np.deg2rad(85)
    angles = angles[(angles >= low_limit) & (angles <= high_limit)]
    spots = np.zeros((len(angles), 2))  # Right side or posterior has n/2 points

    a_se, b_se, n = 2.5, 1.5, 2.2
    na = 2 / n  # Exponent adjustment

    spots[:, 0] = (np.abs(np.cos(angles)) ** na) * a_se * sgn(np.cos(angles))
    spots[:, 1] = (np.abs(np.sin(angles)) ** na) * b_se * sgn(np.sin(angles))

    # Posterior FGs
    gaussian_spaced_points_2 = gauss_points(int(0.4 * n_points))
    angles2 = np.deg2rad(gaussian_spaced_points_2)
    low_limit2 = np.deg2rad(-85)
    high_limit2 = np.deg2rad(+85)
    angles2 = angles2[(angles2 >= low_limit2) & (angles2 <= high_limit2)]
    spots2 = np.zeros((len(angles2), 2))  # Left side or anterior has n/3 points
    spots2[:, 0] = (np.abs(np.cos(angles2)) ** na) * a_se * sgn(np.cos(angles2))
    spots2[:, 1] = (np.abs(np.sin(angles2)) ** na) * b_se * sgn(np.sin(angles2))
    spots = np.concatenate((spots, spots2))

    return spots


def make_fgs_celegans_spindle(motors, a, b, cell, spindle_poles, frame):
    """
    Generate FG spots for C. elegans metaphase/anaphase spindle (posterior
    shift): FGs concentrated more tightly (90deg cutoffs) than in PNC mode.
    """
    n_points = motors

    # Anterior FGs
    gaussian_spaced_points = gauss_points(int(0.4 * n_points))
    angles = np.deg2rad(gaussian_spaced_points) + np.pi
    low_limit = np.deg2rad(90)
    high_limit = 2 * np.pi - np.deg2rad(90)
    angles = angles[(angles >= low_limit) & (angles <= high_limit)]
    spots = np.zeros((len(angles), 2))  # Right side or posterior has n/2 points
    a_se, b_se, n = 2.5, 1.5, 2.2
    na = 2 / n  # Exponent adjustment

    spots[:, 0] = (np.abs(np.cos(angles)) ** na) * a_se * sgn(np.cos(angles))
    spots[:, 1] = (np.abs(np.sin(angles)) ** na) * b_se * sgn(np.sin(angles))
    # Posterior FGs
    gaussian_spaced_points_2 = gauss_points(int(0.6 * n_points))
    angles2 = np.deg2rad(gaussian_spaced_points_2)
    low_limit2 = np.deg2rad(-90)
    high_limit2 = np.deg2rad(+90)
    angles2 = angles2[(angles2 >= low_limit2) & (angles2 <= high_limit2)]
    spots2 = np.zeros((len(angles2), 2))  # Left side or anterior has n/3 points
    spots2[:, 0] = (np.abs(np.cos(angles2)) ** na) * a_se * sgn(np.cos(angles2))
    spots2[:, 1] = (np.abs(np.sin(angles2)) ** na) * b_se * sgn(np.sin(angles2))
    spots = np.concatenate((spots, spots2))

    return spots


def make_fgs_zebrafish_junc(motors, a, b, cell, spindle_poles, frame):
    """
    Generate FG spots for zebrafish junction-only ('junc') mode: FGs
    localized to the cell-cell junction segment tracked from movies.
    """
    juncs = transform_junctions()

    dist1 = distance_matrix(np.array([juncs[frame - 1, 1]]), cell)  # find new closest spot
    end_spot_index = np.argmin(dist1[0])

    dist2 = distance_matrix(np.array([juncs[frame - 1, 0]]), cell)  # find new closest spot
    start_spot_index = np.argmin(dist2[0])

    # Slice FGs from cell array based on junction start and end position and motor density
    partial_perimeter = calculate_partial_perimeter(cell, start_spot_index, end_spot_index)

    n_points_junc = int(motors * partial_perimeter)
    spots = slice_array(cell, start_spot_index, end_spot_index, int(n_points_junc))

    return spots


def make_fgs_zebrafish_uf(motors, a, b, cell, spindle_poles, frame):
    """
    Generate FG spots for zebrafish uniform cortical ('ever'/UF) mode: FGs
    spread uniformly across the whole tracked cell perimeter.
    """
    n_points = int(motors * calculate_perimeter(cell))
    spots = cell[1::int(len(cell) / n_points)]
    return spots


# ============================================================
# Cell boundary / initial spindle placement — cell-type specific
# ============================================================

def make_cell_follicle(params):
    """
    Returns an array (N,2) of x,y coordinates of vertices of a polygon representing the cell
    (follicular epithelium / neuroblast: a simple ellipse).
    """
    a = params[0][0]
    r = params[0][2]  # spindle length=2*r
    FG_density = int(params[1])
    n_astro = int(params[2])
    b = params[0][1]
    spindle_angle = params[3]  # np.random.uniform(0,np.pi)
    theta = np.linspace(0, 2 * np.pi, config.number_of_sides + 1)[:-1].copy()

    cell = np.zeros((config.number_of_sides, 2))
    cell[:, 0] = a * np.cos(theta)
    cell[:, 1] = b * np.sin(theta)

    # Making the spindle
    spindle_poles = np.zeros((2, 2))
    displacement = 1 * np.array([0, 0])

    spindle_poles[0] = [r * np.cos(spindle_angle) + displacement[0], r * np.sin(spindle_angle) + displacement[1]]
    spindle_poles[1] = [r * np.cos(spindle_angle + np.pi) + displacement[0],
                        r * np.sin(spindle_angle + np.pi) + displacement[1]]

    # Sometimes the requested spindle angle doesnt fit so adjustment to the angle is made
    i = 0
    while (point_in_polygon(spindle_poles[0], cell) == False or point_in_polygon(spindle_poles[1], cell) == False):
        spindle_angle = spindle_angle + np.pi / 12
        spindle_poles[0] = [r * np.cos(spindle_angle) + displacement[0], r * np.sin(spindle_angle) + displacement[1]]
        spindle_poles[1] = [r * np.cos(spindle_angle + np.pi) + displacement[0],
                            r * np.sin(spindle_angle + np.pi) + displacement[1]]
        i = i + 1
        if (i > 12):
            break

    # Making motors
    if (NB is not None and NB > 0):
        spots = make_fgs_neuroblast(FG_density, a, b, cell, spindle_poles, 1)
    else:
        spots = make_fgs_follicle(FG_density, a, b, cell, spindle_poles, 1)

    # Making astral microtubules
    astral_MTs, astral_angles, state, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2 = make_astral_MTs(
        params, cell, spindle_poles, spindle_angle, spots)

    return cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, spots, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2


def make_cell_celegans(params):
    """
    Returns an array (N,2) of x,y coordinates of vertices of a polygon representing the cell
    (C. elegans: a superellipse). PNC and spindle modes only differ in the initial
    spindle-pole displacement used for the placement search.
    """
    a = params[0][0]
    r = params[0][2]  # spindle length=2*r
    FG_density = int(params[1])
    n_astro = int(params[2])
    b = params[0][1]
    spindle_angle = params[3]  # np.random.uniform(0,np.pi)
    theta = np.linspace(0, 2 * np.pi, config.number_of_sides + 1)[:-1].copy()

    a, b, n = 2.5, 1.5, 2.2
    na = 2 / n  # Exponent adjustment

    cell = np.zeros((config.number_of_sides, 2))
    cell[:, 0] = (np.abs(np.cos(theta)) ** na) * a * sgn(np.cos(theta))
    cell[:, 1] = (np.abs(np.sin(theta)) ** na) * b * sgn(np.sin(theta))

    # Making the spindle
    spindle_poles = np.zeros((2, 2))

    if model == 'PNC':
        displacement = np.array([1., 0])
    else:
        displacement = np.array([0.0, 0])

    spindle_poles[0] = [r * np.cos(spindle_angle) + displacement[0], r * np.sin(spindle_angle) + displacement[1]]
    spindle_poles[1] = [r * np.cos(spindle_angle + np.pi) + displacement[0],
                        r * np.sin(spindle_angle + np.pi) + displacement[1]]

    # Sometimes the requested spindle angle doesnt fit so adjustment to the angle is made
    i = 0
    while (point_in_polygon(spindle_poles[0], cell) == False or point_in_polygon(spindle_poles[1], cell) == False):
        spindle_angle = spindle_angle + np.pi / 12
        spindle_poles[0] = [r * np.cos(spindle_angle) + displacement[0], r * np.sin(spindle_angle) + displacement[1]]
        spindle_poles[1] = [r * np.cos(spindle_angle + np.pi) + displacement[0],
                            r * np.sin(spindle_angle + np.pi) + displacement[1]]
        i = i + 1
        if (i > 12):
            break

    # Making motors
    if model == 'PNC':
        spots = make_fgs_celegans_pnc(FG_density, a, b, cell, spindle_poles, 1)
    else:
        spots = make_fgs_celegans_spindle(FG_density, a, b, cell, spindle_poles, 1)

    # Making astral microtubules
    astral_MTs, astral_angles, state, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2 = make_astral_MTs(
        params, cell, spindle_poles, spindle_angle, spots)

    return cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, spots, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2


def make_cell_zebrafish(params):
    """
    Returns an array (N,2) of x,y coordinates of vertices of a polygon representing the cell
    (zebrafish: real cell shape tracked from live movies).
    """
    a = params[0][0]
    r = params[0][2]  # spindle length=2*r
    FG_density = int(params[1])
    n_astro = int(params[2])
    b = params[0][1]
    spindle_angle = params[3]  # np.random.uniform(0,np.pi)

    cell = get_real_cell()

    # Making the spindle
    spindle_poles = np.zeros((2, 2))
    displacement = np.array([0, 0])

    spindle_poles[0] = [r * np.cos(spindle_angle) + displacement[0], r * np.sin(spindle_angle) + displacement[1]]
    spindle_poles[1] = [r * np.cos(spindle_angle + np.pi) + displacement[0],
                        r * np.sin(spindle_angle + np.pi) + displacement[1]]

    # Sometimes the requested spindle angle doesnt fit so adjustment to the angle is made
    i = 0
    while (point_in_polygon(spindle_poles[0], cell) == False or point_in_polygon(spindle_poles[1], cell) == False):
        spindle_angle = spindle_angle + np.pi / 12
        spindle_poles[0] = [r * np.cos(spindle_angle) + displacement[0], r * np.sin(spindle_angle) + displacement[1]]
        spindle_poles[1] = [r * np.cos(spindle_angle + np.pi) + displacement[0],
                            r * np.sin(spindle_angle + np.pi) + displacement[1]]
        i = i + 1
        if (i > 12):
            break

    # Making motors
    if model == 'ever':
        spots = make_fgs_zebrafish_uf(FG_density, a, b, cell, spindle_poles, 1)
    elif model == 'junc':
        spots = make_fgs_zebrafish_junc(FG_density, a, b, cell, spindle_poles, 1)

    # Making astral microtubules
    astral_MTs, astral_angles, state, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2 = make_astral_MTs(
        params, cell, spindle_poles, spindle_angle, spots)

    return cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, spots, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2


# ============================================================
# Run-name parsing — cell-type specific (each script's run-name convention
# encodes a slightly different set of parameters)
# ============================================================

def extract_parameters_follicle(input_string):
    motor_density_match = re.search(r'MUD_(\d+)', input_string)
    astral_MTs_match = re.search(r'MT_(\d+)', input_string)
    push_match = re.search(r'push_([\d.]+)', input_string)
    pull_match = re.search(r'pull_([\d.]+)', input_string)
    ts_match = re.search(r'ts_([\d.]+)', input_string)
    cell_match = re.search(r'cell_([A-Za-z]+)', input_string)
    length_MTs_match = re.search(r'AL_([A-Za-z]+)', input_string)
    state_MTs_match = re.search(r'state_([A-Za-z]+)', input_string)
    angle_match = re.search(r'angle_([\d.]+)', input_string)
    SL_match = re.search(r'SL_([\d.]+)', input_string)
    AL_match = re.search(r'AL_([\d.]+)', input_string)
    mcd_match = re.search(r'mcd_([\d.]+)', input_string)
    # Check if in input string 'FE' is followed by 'NB' like in 'FE_NB'
    NB_match = re.search(r'FE_NB', input_string)
    RC_match = re.search(r'RC_rate_([\d.]+)', input_string)
    if cell_match:
        cell_type = cell_match.group(1)  # Extract the matched string ("FE")
        print("Cell type found:", cell_type)  # Output: Cell type found: FE
    else:
        print("No cell type found.")

    if length_MTs_match:
        length_MTs = length_MTs_match.group(1)
    else:
        print("No length_MTs found.")
        length_MTs = "default_length"
    if state_MTs_match:
        state_MTs = state_MTs_match.group(1)
    else:
        print("No state_MTs found.")
        state_MTs = "default_state"
    # Extract values or set to None if not found
    motor_density = int(motor_density_match.group(1)) if motor_density_match else None
    astral_MTs = int(astral_MTs_match.group(1)) if astral_MTs_match else None
    push = float(push_match.group(1)) if push_match else None
    pull = float(pull_match.group(1)) if pull_match else None
    angle = float(angle_match.group(1)) if angle_match else None
    SL = float(SL_match.group(1)) if SL_match else None
    AL = float(AL_match.group(1)) if AL_match else None
    RC_rate = float(RC_match.group(1)) if RC_match else 1.0
    # if NB_match true set NB to 1 else 0
    NB = 1 if NB_match else 0
    mcd = float(mcd_match.group(1)) if mcd_match else None
    time_step = float(ts_match.group(1)) if ts_match else None

    # Return as a dictionary
    return motor_density, astral_MTs, push, pull, SL, time_step, cell_type, angle, NB, RC_rate


def extract_parameters_zebrafish(input_string):
    motor_density_match = re.search(r'MUD_(\d+)', input_string)
    astral_MTs_match = re.search(r'MT_(\d+)', input_string)
    push_match = re.search(r'push_([\d.]+)', input_string)
    pull_match = re.search(r'pull_([\d.]+)', input_string)
    ts_match = re.search(r'ts_([\d.]+)', input_string)
    cell_match = re.search(r'cell_([A-Za-z]+)', input_string)
    length_MTs_match = re.search(r'AL_([A-Za-z]+)', input_string)
    state_MTs_match = re.search(r'state_([A-Za-z]+)', input_string)
    angle_match = re.search(r'angle_([\d.]+)', input_string)
    SL_match = re.search(r'SL_([\d.]+)', input_string)
    AL_match = re.search(r'AL_([\d.]+)', input_string)
    mcd_match = re.search(r'mcd_([\d.]+)', input_string)
    junc_match = re.search(r'junc', input_string)
    ever_match = re.search(r'ever', input_string)

    # Add new regex for GS_rate which is a float number
    GS_rate_match = re.search(r'GS_rate_([\d.]+)', input_string)
    if cell_match:
        cell_type = cell_match.group(1)  # Extract the matched string ("FE")
        print("Cell type found:", cell_type)  # Output: Cell type found: FE
    else:
        print("No cell type found.")

    if length_MTs_match:
        length_MTs = length_MTs_match.group(1)
    else:
        print("No length_MTs found.")
        length_MTs = "default_length"
    if state_MTs_match:
        state_MTs = state_MTs_match.group(1)
    else:
        print("No state_MTs found.")
        state_MTs = "default_state"
    # Extract values or set to None if not found
    motor_density = int(motor_density_match.group(1)) if motor_density_match else None
    astral_MTs = int(astral_MTs_match.group(1)) if astral_MTs_match else None
    push = float(push_match.group(1)) if push_match else None
    pull = float(pull_match.group(1)) if pull_match else None
    angle = float(angle_match.group(1)) if angle_match else None
    SL = float(SL_match.group(1)) if SL_match else None
    AL = float(AL_match.group(1)) if AL_match else 1
    # set GS_rate to 1 if not found, otherwise extract the value
    GS_rate = float(GS_rate_match.group(1)) if GS_rate_match else 1.0
    # set model to 'ever' if ever is found, if found junc set it to junc
    if ever_match:
        model = 'ever'
    if junc_match:
        model = 'junc'
    mcd = float(mcd_match.group(1)) if mcd_match else None
    time_step = float(ts_match.group(1)) if ts_match else None

    # Return as a dictionary
    return motor_density, astral_MTs, push, pull, SL, time_step, cell_type, angle, model, GS_rate, AL


def extract_parameters_celegans(input_string):
    motor_density_match = re.search(r'MUD_(\d+)', input_string)
    astral_MTs_match = re.search(r'MT_(\d+)', input_string)
    push_match = re.search(r'push_([\d.]+)', input_string)
    pull_match = re.search(r'pull_([\d.]+)', input_string)
    ts_match = re.search(r'ts_([\d.]+)', input_string)
    cell_match = re.search(r'cell_([A-Za-z]+)', input_string)
    length_MTs_match = re.search(r'AL_([A-Za-z]+)', input_string)
    state_MTs_match = re.search(r'state_([A-Za-z]+)', input_string)
    angle_match = re.search(r'angle_([\d.]+)', input_string)
    SL_match = re.search(r'SL_([\d.]+)', input_string)
    AL_match = re.search(r'AL_([\d.]+)', input_string)
    mcd_match = re.search(r'mcd_([\d.]+)', input_string)
    if cell_match:
        cell_type = cell_match.group(1)  # Extract the matched string ("FE")
        print("Cell type found:", cell_type)  # Output: Cell type found: FE
    else:
        print("No cell type found.")

    if length_MTs_match:
        length_MTs = length_MTs_match.group(1)
    else:
        print("No length_MTs found.")
        length_MTs = "default_length"

    if state_MTs_match:
        state_MTs = state_MTs_match.group(1)
    else:
        print("No state_MTs found.")
        state_MTs = "default_state"
    # Extract values or set to None if not found
    if "_PNC_" in input_string:
        model = "PNC"
        angle = 90
    elif "_position_" in input_string:
        model = "spindle"
        angle = 0
    motor_density = int(motor_density_match.group(1)) if motor_density_match else None
    astral_MTs = int(astral_MTs_match.group(1)) if astral_MTs_match else None
    push = float(push_match.group(1)) if push_match else None
    pull = float(pull_match.group(1)) if pull_match else None

    SL = float(SL_match.group(1)) if SL_match else None
    AL = float(AL_match.group(1)) if AL_match else None
    mcd = float(mcd_match.group(1)) if mcd_match else None
    time_step = float(ts_match.group(1)) if ts_match else None

    # Return as a dictionary
    return motor_density, astral_MTs, push, pull, SL, time_step, cell_type, angle, model


# ============================================================
# Per-frame plotting — cell-type specific (purely cosmetic; does not affect
# simulation results)
# ============================================================

def plot_cell_follicle(cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, push_force, pull_force,
                       spots, i, which_bind, which_push, v_c, n, ratio, params, borders):
    """The tricky part I sometimes forget about is that all the data in the legend such as total force, total pushing force and pulling force are for current iteration i,
    but all the vectors (lines) are shown for future iteration. Basically, vectors predict the next step."""
    # Initialization
    tpoint = i
    a = params[0][0]
    b = params[0][1]
    left_min = min(cell[:, 0]) - 1.5 * a
    right_max = max(cell[:, 0]) + 1.5 * a
    top = max(cell[:, 1]) + 0.5
    bottom = min(cell[:, 1]) - 0.5

    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.set_xlim([left_min, right_max])
    ax.set_ylim([bottom, top])

    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_visible(True)

    # Calculate force vectors for display
    count_pushing = np.sum(which_push == 1, axis=1)
    count_pulling = np.sum(which_bind == 1, axis=1)
    center = np.array(
        [1 / 2 * (spindle_poles[0, 0] + spindle_poles[1, 0]), 1 / 2 * (spindle_poles[0, 1] + spindle_poles[1, 1])])

    # Calculating repel force
    dist_to_cort = distance_matrix(spindle_poles, cell)
    repel_force1 = np.array([0, 0])
    repel_force2 = np.array([0, 0])

    if (np.min(dist_to_cort[0]) <= 1.5 * config.min_cortex_dist):
        repel_vec1 = -spindle_poles[0] / LA.norm(spindle_poles[0])
        repel_force1 = (config.repel / np.min(dist_to_cort[0])) * repel_vec1

    if (np.min(dist_to_cort[1]) <= 1.5 * config.min_cortex_dist):
        repel_vec2 = -spindle_poles[1] / LA.norm(spindle_poles[1])
        repel_force2 = (config.repel / np.min(dist_to_cort[1])) * repel_vec2

    pull_t = np.sum(pull_force[0] + pull_force[1], axis=0)
    push_t = np.sum(push_force[0] + push_force[1], axis=0)

    force_vector_1 = pull_force[0] + push_force[0]
    force_vector_2 = pull_force[1] + push_force[1]
    force1 = np.sum(force_vector_1, axis=0) + repel_force1
    force2 = np.sum(force_vector_2, axis=0) + repel_force2
    pull_1_t = np.sum(pull_force[0], axis=0)
    push_1_t = np.sum(push_force[0], axis=0)
    pull_2_t = np.sum(pull_force[1], axis=0)
    push_2_t = np.sum(push_force[1], axis=0)

    if (LA.norm(force1 + force2) == 0):
        f_net = np.array([0, 0])
    else:
        f_net = (force1 + force2) / LA.norm(force1 + force2)

    if ((LA.norm(pull_t) + LA.norm(push_t)) == 0):
        ratio = 0
    else:
        ratio = 100 * LA.norm(pull_t) / (LA.norm(pull_t) + LA.norm(push_t))

    ax.plot(cell[:, 0], cell[:, 1], color='darkgrey', label='Time = %.1f s' % (n * time_step), linewidth=5, zorder=10)

    ax.scatter([spindle_poles[0, 0], spindle_poles[1, 0]], [spindle_poles[0, 1], spindle_poles[1, 1]], color='yellow',
               s=150, edgecolors='k', zorder=6)

    # CORTEX
    theta = np.linspace(0, 2 * np.pi, config.number_of_sides + 1)[:-1].copy()
    cortex = np.zeros((config.number_of_sides, 2))
    cortex[:, 0] = 1.02 * a * np.cos(theta)
    cortex[:, 1] = 1.02 * b * np.sin(theta)
    ax.plot([cell[0, 0], cell[-1, 0]],
            [cell[0, 1], cell[-1, 1]], color='darkgrey', label=f'Angle ={math.degrees(spindle_angle):.2f}°')

    # MOTORS
    if pull != 0:
        ax.scatter(spots[:, 0], spots[:, 1],
                   label='Number of astro MTs = %.1f,\n Number of FGs = %.1f' % (2 * len(astral_MTs[0]), len(spots)),
                   color='salmon', s=55, edgecolors='dimgrey', marker="8", zorder=50)

    # COLORCODING MTs
    for i in range(len(astral_MTs[0])):
        if (which_bind[0, i] == 1):  # binded
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'tab:red', zorder=5)
        elif (state[0, i] == -1):  # shrinking
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'tab:cyan', zorder=3)
        elif (which_push[0, i] == 1):  # pushing
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'darkgreen', zorder=4)
        else:  # which_push=0, state=1,which_bind=0
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'slateblue', zorder=2)

    for i in range(len(astral_MTs[1])):
        if (which_bind[1, i] == 1):
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'tab:red', zorder=5)
        elif (state[1, i] == -1):
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'tab:cyan', zorder=3)
        elif (which_push[1, i] == 1):
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'darkgreen', zorder=4)
        else:
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'slateblue', zorder=2)

    # Force vectors-------------------------------------------------------------------------------------------
    com = np.array([(spindle_poles[0, 0] + spindle_poles[1, 0]) / 2, (spindle_poles[0, 1] + spindle_poles[1, 1]) / 2])

    # ANNOTATION
    if (config.show_vectors == 1):
        kap = 0.01
        # EACH POLE SEPARATE PULL AND PUSH
        ax.plot([spindle_poles[0, 0], spindle_poles[0, 0] - kap * push_1_t[0]],
                [spindle_poles[0, 1], spindle_poles[0, 1] - kap * push_1_t[1]], 'tab:red', linewidth=4,
                label='F_push 1 = %.2f' % (LA.norm(push_1_t)))
        ax.plot([spindle_poles[0, 0], spindle_poles[0, 0] + kap * pull_1_t[0]],
                [spindle_poles[0, 1], spindle_poles[0, 1] + kap * pull_1_t[1]], 'tab:pink', linewidth=4,
                label='F_pull 1 = %.2f' % (LA.norm(pull_1_t)))
        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] - kap * push_2_t[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] - kap * push_2_t[1]], 'tab:purple', linewidth=4,
                label='F_push 2 = %.2f' % (LA.norm(push_2_t)))
        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] + kap * pull_2_t[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] + kap * pull_2_t[1]], 'tab:blue', linewidth=4,
                label='F_pull 2 = %.2f' % (LA.norm(pull_2_t)))

        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] - kap * repel_force1[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] - kap * repel_force1[1]], 'indigo', linewidth=4,
                label='repel 1 = %.2f' % (LA.norm(repel_force1)))
        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] + kap * repel_force2[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] + kap * repel_force2[1]], 'slategray', linewidth=4,
                label='repel 2 = %.2f' % (LA.norm(repel_force2)))
    else:
        kap = 0
        fig.patch.set_visible(False)
        ax.axis('off')

    # NET FORCE EACH POLE
    ax.plot([spindle_poles[0, 0], spindle_poles[0, 0] + kap * force1[0]],
            [spindle_poles[0, 1], spindle_poles[0, 1] + kap * force1[1]], 'tab:orange', ls='--', linewidth=4,
            label=r'$F_{{\text{{pole1}}}} ={:.3f} \ pN $'.format(LA.norm(force1)), zorder=90)
    ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] + kap * force2[0]],
            [spindle_poles[1, 1], spindle_poles[1, 1] + kap * force2[1]], 'k', ls='--', linewidth=4,
            label=r'$F_{{\text{{pole2}}}} ={:.3f} \ pN $'.format(LA.norm(force2)), zorder=90)

    # NET FORCE
    ax.plot([com[0], com[0] + kap * f_net[0]], [com[1], com[1] + kap * f_net[1]], 'k', linewidth=4,
            label='Total force = %.2f pN' % (LA.norm(force1 + force2)))

    # TEXT
    plt.text(right_max - 0.75, top - 0.25, r'$F_{{\text{{pull}}}} / F_{{\text{{push}}}} = {:.3f}$'.format(ratio),
             fontsize=10)
    variable1 = np.sum(count_pulling)
    variable2 = np.sum(count_pushing)
    plt.text(right_max - 0.75, top - 0.35,
             r'$N_{{\text{{pull}}}} = {}, N_{{\text{{push}}}} = {}$'.format(variable1, variable2), fontsize=10)
    rect = patches.Rectangle((right_max - 0.6, top - 0.85), 0.15, 0.15, facecolor=create_color_gradient(ratio / 100))
    ax.add_patch(rect)

    # Parameters for the free body diagram
    if (config.show_vectors == 1):
        rod_length = 0.6  # Length of the rod
        rod_angle = np.rad2deg(spindle_angle)  # Angle of the rod in degrees
        sphere_radius = 0.1  # Radius of the spheres at the ends of the rod
        rod_angle_rad = np.radians(rod_angle)
        center_x, center_y = 1.8, -0.2  # Center of the rod
        end1_x = center_x + (rod_length / 2) * np.cos(rod_angle_rad)
        end1_y = center_y + (rod_length / 2) * np.sin(rod_angle_rad)
        end2_x = center_x - (rod_length / 2) * np.cos(rod_angle_rad)
        end2_y = center_y - (rod_length / 2) * np.sin(rod_angle_rad)

        rod = patches.FancyArrowPatch((end1_x, end1_y), (end2_x, end2_y),
                                      arrowstyle='-', color='black', lw=2)
        ax.add_patch(rod)

        sphere1 = patches.Circle((end1_x, end1_y), sphere_radius, color='gold')
        sphere2 = patches.Circle((end2_x, end2_y), sphere_radius, color='gold')
        ax.add_patch(sphere1)
        ax.add_patch(sphere2)

        forces1 = np.array([push_1_t, pull_1_t, repel_force1, force1])
        forces2 = np.array([push_2_t, pull_2_t, repel_force2, force2])

        colors1 = ['tab:red', 'tab:pink', 'green', 'black']
        colors2 = ['tab:purple', 'blue', 'green', 'black']

        labels1 = ["Push1", "Pull1", "Repel1", "Net1"]
        labels2 = ["Push2", "Pull2", "Repel2", "Net2"]

        draw_vectors(ax, (end1_x, end1_y), forces1, labels1, colors1)
        draw_vectors(ax, (end2_x, end2_y), forces2, labels2, colors2)

    # CHROMOSOMES
    chrom_angle = np.pi / 2 + spindle_angle
    n_chrom = 4
    dd = 0.8 * w
    xs = np.linspace(com[0] - dd * np.cos(chrom_angle), com[0] + dd * np.cos(chrom_angle), n_chrom)
    ys = np.linspace(com[1] - dd * np.sin(chrom_angle), com[1] + dd * np.sin(chrom_angle), n_chrom)

    for i in range(n_chrom):
        line1, line2 = chrom(xs[i], ys[i], chrom_angle)
        plt.plot([spindle_poles[0, 0], xs[i]], [spindle_poles[0, 1], ys[i]], color='g', linewidth=10)
        plt.plot([spindle_poles[1, 0], xs[i]], [spindle_poles[1, 1], ys[i]], color='g', linewidth=10)
        plt.plot(line1[0], line1[1], color='dodgerblue', linewidth=14)
        plt.plot(line2[0], line2[1], color='dodgerblue', linewidth=14)

    ax.legend(ncol=1, loc='upper left', facecolor='white', framealpha=1)

    # SAVING
    path = params[5]
    name = params[4]
    plt.savefig(path + '/' + name + '_' + str(1 + n) + '.pdf', bbox_inches='tight', pad_inches=1)
    plt.close(fig)


# ============================================================
# Per-frame plotting — C. elegans
# ============================================================
def plot_cell_celegans(cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, push_force, pull_force,
                       spots, i, which_bind, which_push, v_c, n, ratio, params, borders):
    """The tricky part I sometimes forget about is that all the data in the legend such as total force, total pushing force and pulling force are for current iteration i,
    but all the vectors (lines) are shown for future iteration. Basically, vectors predict the next step."""
    tpoint = i
    a = params[0][0]
    b = params[0][1]
    r = LA.norm(spindle_poles[0] - spindle_poles[1]) / 2  # spindle length=2*r
    left_min = min(cell[:, 0]) - .5 * a
    right_max = max(cell[:, 0]) + .5 * a
    top = max(cell[:, 1]) + 1.
    bottom = min(cell[:, 1]) - 1.

    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.set_xlim([left_min, right_max])
    ax.set_ylim([bottom, top])

    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_visible(True)

    count_pushing = np.sum(which_push == 1, axis=1)
    count_pulling = np.sum(which_bind == 1, axis=1)
    center = np.array(
        [1 / 2 * (spindle_poles[0, 0] + spindle_poles[1, 0]), 1 / 2 * (spindle_poles[0, 1] + spindle_poles[1, 1])])

    dist_to_cort = distance_matrix(spindle_poles, cell)
    repel_force1 = np.array([0, 0])
    repel_force2 = np.array([0, 0])

    if (np.min(dist_to_cort[0]) <= 1.5 * config.min_cortex_dist):
        repel_vec1 = -spindle_poles[0] / LA.norm(spindle_poles[0])
        repel_force1 = (config.repel / np.min(dist_to_cort[0])) * repel_vec1

    if (np.min(dist_to_cort[1]) <= 1.5 * config.min_cortex_dist):
        repel_vec2 = -spindle_poles[1] / LA.norm(spindle_poles[1])
        repel_force2 = (config.repel / np.min(dist_to_cort[1])) * repel_vec2

    pull_t = np.sum(pull_force[0] + pull_force[1], axis=0)
    push_t = np.sum(push_force[0] + push_force[1], axis=0)

    force_vector_1 = pull_force[0] + push_force[0]
    force_vector_2 = pull_force[1] + push_force[1]
    force1 = np.sum(force_vector_1, axis=0) + repel_force1
    force2 = np.sum(force_vector_2, axis=0) + repel_force2
    pull_1_t = np.sum(pull_force[0], axis=0)
    push_1_t = np.sum(push_force[0], axis=0)
    pull_2_t = np.sum(pull_force[1], axis=0)
    push_2_t = np.sum(push_force[1], axis=0)

    if (LA.norm(force1 + force2) == 0):
        f_net = np.array([0, 0])
    else:
        f_net = (force1 + force2) / LA.norm(force1 + force2)

    if ((LA.norm(pull_t) + LA.norm(push_t)) == 0):
        ratio = 0
    else:
        ratio = 100 * LA.norm(pull_t) / (LA.norm(pull_t) + LA.norm(push_t))

    ax.plot(cell[:, 0], cell[:, 1], color='dimgrey', label='Time = %.1f s' % (n * time_step), linewidth=5, zorder=10)

    ax.scatter([spindle_poles[0, 0], spindle_poles[1, 0]], [spindle_poles[0, 1], spindle_poles[1, 1]], color='yellow',
               s=150, edgecolors='k', zorder=6, label=f'Spindle length = {20 * r:.5f} μm')
    if model == 'PNC':
        ax.plot([spindle_poles[0, 0], spindle_poles[1, 0]], [spindle_poles[0, 1], spindle_poles[1, 1]], color='grey',
                linewidth=7, zorder=5)

    # CORTEX
    theta = np.linspace(0, 2 * np.pi, config.number_of_sides + 1)[:-1].copy()
    cortex = np.zeros((config.number_of_sides, 2))
    cortex[:, 0] = 1.02 * a * np.cos(theta)
    cortex[:, 1] = 1.02 * b * np.sin(theta)
    ax.plot([cell[0, 0], cell[-1, 0]],
            [cell[0, 1], cell[-1, 1]], color='darkgrey', label=f'Angle ={math.degrees(spindle_angle):.2f}°')

    # MOTORS
    ax.scatter(spots[:, 0], spots[:, 1],
               label='Number of astro MTs = %.1f,\n Number of FGs = %.1f' % (2 * len(astral_MTs[0]), len(spots)),
               color='salmon', s=55, edgecolors='dimgrey', marker="8", zorder=50)

    # COLORCODING MTs
    for i in range(len(astral_MTs[0])):
        if (which_bind[0, i] == 1):  # binded
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'red', zorder=5)
        elif (state[0, i] == -1):  # shrinking
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'tab:cyan', zorder=3)
        elif (which_push[0, i] == 1):  # pushing
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'darkgreen', zorder=4)
        else:  # which_push=0, state=1,which_bind=0
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'slateblue', zorder=2)

    for i in range(len(astral_MTs[1])):
        if (which_bind[1, i] == 1):
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'tab:red', zorder=5)
        elif (state[1, i] == -1):
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'tab:cyan', zorder=3)
        elif (which_push[1, i] == 1):
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'green', zorder=4)
        else:
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'slateblue', zorder=2)

    # Force vectors-------------------------------------------------------------------------------------------
    spindle_envelope = np.zeros((359, 2))
    theta = np.linspace(0, 2 * np.pi, 360)[:-1].copy()
    com = np.array([(spindle_poles[0, 0] + spindle_poles[1, 0]) / 2, (spindle_poles[0, 1] + spindle_poles[1, 1]) / 2])
    spindle_envelope[:, 0] = r * np.cos(theta) + com[0]
    spindle_envelope[:, 1] = w * np.sin(theta) + com[-1]
    spindle_envelope = rotate_points(spindle_envelope, spindle_angle)
    if (model == 'PNC'):
        ax.plot(spindle_envelope[:, 0], spindle_envelope[:, 1], linewidth=7, color='grey', zorder=1)
        plt.fill(spindle_envelope[:, 0], spindle_envelope[:, 1], color='lightsteelblue', edgecolor='darkgrey',
                 linewidth=2)

    # MEAN MTs LENGTH
    leng = np.zeros((2, len(state[0])))
    for i in range(2):
        for j in range(len(astral_MTs[0])):
            leng[i, j] = LA.norm(
                [astral_MTs[i, j, -1, 0] - spindle_poles[i, 0], astral_MTs[i, j, -1, 1] - spindle_poles[i, 1]])
    mean_length = np.average(leng)

    # ANNOTATION
    if (config.show_vectors == 1):
        kap = 0.01
        ax.plot([spindle_poles[0, 0], spindle_poles[0, 0] - kap * push_1_t[0]],
                [spindle_poles[0, 1], spindle_poles[0, 1] - kap * push_1_t[1]], 'tab:red', linewidth=4,
                label='F_push 1 = %.2f' % (LA.norm(push_1_t)))
        ax.plot([spindle_poles[0, 0], spindle_poles[0, 0] + kap * pull_1_t[0]],
                [spindle_poles[0, 1], spindle_poles[0, 1] + kap * pull_1_t[1]], 'tab:pink', linewidth=4,
                label='F_pull 1 = %.2f' % (LA.norm(pull_1_t)))
        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] - kap * push_2_t[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] - kap * push_2_t[1]], 'tab:purple', linewidth=4,
                label='F_push 2 = %.2f' % (LA.norm(push_2_t)))
        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] + kap * pull_2_t[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] + kap * pull_2_t[1]], 'tab:blue', linewidth=4,
                label='F_pull 2 = %.2f' % (LA.norm(pull_2_t)))

        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] - kap * repel_force1[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] - kap * repel_force1[1]], 'indigo', linewidth=4,
                label='repel 1 = %.2f' % (LA.norm(repel_force1)))
        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] + kap * repel_force2[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] + kap * repel_force2[1]], 'slategray', linewidth=4,
                label='repel 2 = %.2f' % (LA.norm(repel_force2)))
    else:
        kap = 0
        fig.patch.set_visible(False)
        ax.axis('off')

    # TEXT
    plt.text(right_max - 2.1, top - 0.3, 'Pull/Push = %.2f' % (ratio), fontsize=12)
    plt.text(right_max - 2.1, top - 0.4, f'Pull,Push = {np.sum(count_pulling), np.sum(count_pushing)}', fontsize=12)
    rect = patches.Rectangle((right_max - 1.2, top - 1.1), 0.3, 0.3, facecolor=create_color_gradient(ratio / 100))
    ax.add_patch(rect)

    # Parameters for the free body diagram
    if (config.show_vectors == 3):
        rod_length = 0.6  # Length of the rod
        rod_angle = np.rad2deg(spindle_angle)  # Angle of the rod in degrees
        sphere_radius = 0.1  # Radius of the spheres at the ends of the rod
        rod_angle_rad = np.radians(rod_angle)
        center_x, center_y = 1.8, -0.2  # Center of the rod
        end1_x = center_x + (rod_length / 2) * np.cos(rod_angle_rad)
        end1_y = center_y + (rod_length / 2) * np.sin(rod_angle_rad)
        end2_x = center_x - (rod_length / 2) * np.cos(rod_angle_rad)
        end2_y = center_y - (rod_length / 2) * np.sin(rod_angle_rad)

        rod = patches.FancyArrowPatch((end1_x, end1_y), (end2_x, end2_y),
                                      arrowstyle='-', color='black', lw=2)
        ax.add_patch(rod)

        sphere1 = patches.Circle((end1_x, end1_y), sphere_radius, color='gold')
        sphere2 = patches.Circle((end2_x, end2_y), sphere_radius, color='gold')
        ax.add_patch(sphere1)
        ax.add_patch(sphere2)

        forces1 = np.array([push_1_t, pull_1_t, repel_force1, force1])
        forces2 = np.array([push_2_t, pull_2_t, repel_force2, force2])

        colors1 = ['tab:red', 'tab:pink', 'green', 'black']
        colors2 = ['tab:purple', 'blue', 'green', 'black']

        labels1 = ["Push1", "Pull1", "Repel1", "Net1"]
        labels2 = ["Push2", "Pull2", "Repel2", "Net2"]

        draw_vectors(ax, (end1_x, end1_y), forces1, labels1, colors1)
        draw_vectors(ax, (end2_x, end2_y), forces2, labels2, colors2)

    # CHROMOSOMES / SEVERED-SPINDLE MARKER
    if (model == 'spindle' and tpoint * time_step < severance_time):
        chrom_angle = np.pi / 2 + spindle_angle
        n_chrom = 4
        dd = 0.5 * w
        xs = np.linspace(com[0] - dd * np.cos(chrom_angle), com[0] + dd * np.cos(chrom_angle), n_chrom)
        ys = np.linspace(com[1] - dd * np.sin(chrom_angle), com[1] + dd * np.sin(chrom_angle), n_chrom)

        for i in range(n_chrom):
            line1, line2 = chrom(xs[i], ys[i], chrom_angle)
            plt.plot([spindle_poles[0, 0], xs[i]], [spindle_poles[0, 1], ys[i]], color='g', linewidth=6)
            plt.plot([spindle_poles[1, 0], xs[i]], [spindle_poles[1, 1], ys[i]], color='g', linewidth=6)
            plt.plot(line1[0], line1[1], color='dodgerblue', linewidth=10)
            plt.plot(line2[0], line2[1], color='dodgerblue', linewidth=10)
    elif (severance == 1 and tpoint * time_step >= severance_time):
        plt.plot([0 - 1.2 * w * np.cos(np.deg2rad(init_sp_angle) + np.pi / 2),
                  0 + 1.2 * w * np.cos(np.deg2rad(init_sp_angle) + np.pi / 2)],
                 [0 - 1.2 * w * np.sin(np.deg2rad(init_sp_angle) + np.pi / 2),
                  0 + 1.2 * w * np.sin(np.deg2rad(init_sp_angle) + np.pi / 2)], color='r', linewidth=15)

    ax.legend(ncol=1, loc='upper left', facecolor='white', framealpha=1)

    # SAVING
    path = params[5]
    name = params[4]
    plt.savefig(path + '/' + name + '_' + str(1 + n) + '.pdf', bbox_inches='tight', pad_inches=1)
    plt.close(fig)


# ============================================================
# Per-frame plotting — zebrafish
# ============================================================
def plot_cell_zebrafish(cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, push_force, pull_force,
                        spots, i, which_bind, which_push, v_c, n, ratio, params, borders):
    """The tricky part I sometimes forget about is that all the data in the legend such as total force, total pushing force and pulling force are for current iteration i,
    but all the vectors (lines) are shown for future iteration. Basically, vectors predict the next step."""
    # Initialization
    tpoint = i
    a = params[0][0]
    b = params[0][1]
    left_min = borders[0]
    right_max = borders[1]
    top = borders[2]
    bottom = borders[3]

    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.set_xlim([left_min, right_max])
    ax.set_ylim([bottom, top])

    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_visible(True)

    # Calculate force vectors for display
    count_pushing = np.sum(which_push == 1, axis=1)
    count_pulling = np.sum(which_bind == 1, axis=1)
    center = np.array(
        [1 / 2 * (spindle_poles[0, 0] + spindle_poles[1, 0]), 1 / 2 * (spindle_poles[0, 1] + spindle_poles[1, 1])])

    # Calculating repel force
    dist_to_cort = distance_matrix(spindle_poles, cell)
    repel_force1 = np.array([0, 0])
    repel_force2 = np.array([0, 0])

    if (np.min(dist_to_cort[0]) <= 1.5 * config.min_cortex_dist):
        repel_vec1 = -spindle_poles[0] / LA.norm(spindle_poles[0])
        repel_force1 = (config.repel / np.min(dist_to_cort[0])) * repel_vec1

    if (np.min(dist_to_cort[1]) <= 1.5 * config.min_cortex_dist):
        repel_vec2 = -spindle_poles[1] / LA.norm(spindle_poles[1])
        repel_force2 = (config.repel / np.min(dist_to_cort[1])) * repel_vec2

    pull_t = np.sum(pull_force[0] + pull_force[1], axis=0)
    push_t = np.sum(push_force[0] + push_force[1], axis=0)

    force_vector_1 = pull_force[0] + push_force[0]
    force_vector_2 = pull_force[1] + push_force[1]
    force1 = np.sum(force_vector_1, axis=0) + repel_force1
    force2 = np.sum(force_vector_2, axis=0) + repel_force2
    pull_1_t = np.sum(pull_force[0], axis=0)
    push_1_t = np.sum(push_force[0], axis=0)
    pull_2_t = np.sum(pull_force[1], axis=0)
    push_2_t = np.sum(push_force[1], axis=0)

    if (LA.norm(force1 + force2) == 0):
        f_net = np.array([0, 0])
    else:
        f_net = (force1 + force2) / LA.norm(force1 + force2)

    if ((LA.norm(pull_t) + LA.norm(push_t)) == 0):
        ratio = 0
    else:
        ratio = 100 * LA.norm(pull_t) / (LA.norm(pull_t) + LA.norm(push_t))

    ax.plot(cell[:, 0], cell[:, 1], color='dimgrey', label='Time = %.1f s' % (n * time_step), linewidth=5, zorder=10)
    ax.legend(prop={'size': 14}, facecolor='white', framealpha=1)

    ax.scatter([spindle_poles[0, 0], spindle_poles[1, 0]], [spindle_poles[0, 1], spindle_poles[1, 1]], color='yellow',
               s=150, edgecolors='k', zorder=6)

    # CORTEX
    theta = np.linspace(0, 2 * np.pi, config.number_of_sides + 1)[:-1].copy()
    cortex = np.zeros((config.number_of_sides, 2))
    cortex[:, 0] = 1.02 * a * np.cos(theta)
    cortex[:, 1] = 1.02 * b * np.sin(theta)
    ax.plot([cell[0, 0], cell[-1, 0]],
            [cell[0, 1], cell[-1, 1]], color='dimgrey',
            label=f'Spindle angle ={math.degrees(abs(long_axis - spindle_angle)):.2f}°')

    # MOTORS
    if pull != 0:
        ax.scatter(spots[:, 0], spots[:, 1],
                   label='Number of astro MTs = %.1f,\n FG density = %.2f' % (2 * len(astral_MTs[0]), params[1]),
                   color='salmon', s=55, edgecolors='dimgrey', marker="8", zorder=50)

    # COLORCODING MTs
    for i in range(len(astral_MTs[0])):
        if (which_bind[0, i] == 1):  # binded
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'tab:red', zorder=5)
        elif (state[0, i] == -1):  # shrinking
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'tab:cyan', zorder=3)
        elif (which_push[0, i] == 1):  # pushing
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'darkgreen', zorder=4)
        else:  # which_push=0, state=1,which_bind=0
            ax.plot(astral_MTs[0, i, :, 0], astral_MTs[0, i, :, 1], 'slateblue', zorder=2)

    for i in range(len(astral_MTs[1])):
        if (which_bind[1, i] == 1):
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'tab:red', zorder=5)
        elif (state[1, i] == -1):
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'tab:cyan', zorder=3)
        elif (which_push[1, i] == 1):
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'darkgreen', zorder=4)
        else:
            ax.plot(astral_MTs[1, i, :, 0], astral_MTs[1, i, :, 1], 'slateblue', zorder=2)

    # Force vectors-------------------------------------------------------------------------------------------
    com = np.array([(spindle_poles[0, 0] + spindle_poles[1, 0]) / 2, (spindle_poles[0, 1] + spindle_poles[1, 1]) / 2])
    spindle_envelope = generate_spindle(spindle_poles, spindle_angle, r, w)

    # ANNOTATION
    if (config.show_vectors == 100):
        kap = 0.01
        # EACH POLE SEPARATE PULL AND PUSH
        ax.plot([spindle_poles[0, 0], spindle_poles[0, 0] - kap * push_1_t[0]],
                [spindle_poles[0, 1], spindle_poles[0, 1] - kap * push_1_t[1]], 'tab:red', linewidth=4,
                label='F_push 1 = %.2f' % (LA.norm(push_1_t)))
        ax.plot([spindle_poles[0, 0], spindle_poles[0, 0] + kap * pull_1_t[0]],
                [spindle_poles[0, 1], spindle_poles[0, 1] + kap * pull_1_t[1]], 'tab:pink', linewidth=4,
                label='F_pull 1 = %.2f' % (LA.norm(pull_1_t)))
        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] - kap * push_2_t[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] - kap * push_2_t[1]], 'tab:purple', linewidth=4,
                label='F_push 2 = %.2f' % (LA.norm(push_2_t)))
        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] + kap * pull_2_t[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] + kap * pull_2_t[1]], 'tab:blue', linewidth=4,
                label='F_pull 2 = %.2f' % (LA.norm(pull_2_t)))

        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] - kap * repel_force1[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] - kap * repel_force1[1]], 'indigo', linewidth=4,
                label='repel 1 = %.2f' % (LA.norm(repel_force1)))
        ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] + kap * repel_force2[0]],
                [spindle_poles[1, 1], spindle_poles[1, 1] + kap * repel_force2[1]], 'slategray', linewidth=4,
                label='repel 2 = %.2f' % (LA.norm(repel_force2)))
    else:
        kap = 0
        fig.patch.set_visible(False)
        ax.axis('off')

    # NET FORCE EACH POLE
    ax.plot([spindle_poles[0, 0], spindle_poles[0, 0] + kap * force1[0]],
            [spindle_poles[0, 1], spindle_poles[0, 1] + kap * force1[1]], 'tab:orange', ls='--', linewidth=4,
            label=r'$F_{{\text{{pole1}}}} ={:.3f} \ pN $'.format(LA.norm(force1)), zorder=90)
    ax.plot([spindle_poles[1, 0], spindle_poles[1, 0] + kap * force2[0]],
            [spindle_poles[1, 1], spindle_poles[1, 1] + kap * force2[1]], 'k', ls='--', linewidth=4,
            label=r'$F_{{\text{{pole2}}}} ={:.3f} \ pN $'.format(LA.norm(force2)), zorder=90)

    # TEXT
    plt.text(right_max - 5.2, top - 2.3, 'Pull/Push = %.2f' % (ratio), fontsize=12)
    plt.text(right_max - 5.2, top - 2.4, f'Pull,Push = {np.sum(count_pulling), np.sum(count_pushing)}', fontsize=12)

    # Parameters for the free body diagram
    if (config.show_vectors == 2):
        rod_length = 0.6  # Length of the rod
        rod_angle = np.rad2deg(spindle_angle)  # Angle of the rod in degrees
        sphere_radius = 0.1  # Radius of the spheres at the ends of the rod
        rod_angle_rad = np.radians(rod_angle)
        center_x, center_y = 1.8, -0.2  # Center of the rod
        end1_x = center_x + (rod_length / 2) * np.cos(rod_angle_rad)
        end1_y = center_y + (rod_length / 2) * np.sin(rod_angle_rad)
        end2_x = center_x - (rod_length / 2) * np.cos(rod_angle_rad)
        end2_y = center_y - (rod_length / 2) * np.sin(rod_angle_rad)

        rod = patches.FancyArrowPatch((end1_x, end1_y), (end2_x, end2_y),
                                      arrowstyle='-', color='black', lw=2)
        ax.add_patch(rod)

        sphere1 = patches.Circle((end1_x, end1_y), sphere_radius, color='gold')
        sphere2 = patches.Circle((end2_x, end2_y), sphere_radius, color='gold')
        ax.add_patch(sphere1)
        ax.add_patch(sphere2)

        forces1 = np.array([push_1_t, pull_1_t, repel_force1, force1])
        forces2 = np.array([push_2_t, pull_2_t, repel_force2, force2])

        colors1 = ['tab:red', 'tab:pink', 'green', 'black']
        colors2 = ['tab:purple', 'blue', 'green', 'black']

        labels1 = ["Push1", "Pull1", "Repel1", "Net1"]
        labels2 = ["Push2", "Pull2", "Repel2", "Net2"]

        draw_vectors(ax, (end1_x, end1_y), forces1, labels1, colors1)
        draw_vectors(ax, (end2_x, end2_y), forces2, labels2, colors2)

    # CHROMOSOMES
    if (cell_type != 'celegans'):
        chrom_angle = np.pi / 2 + spindle_angle
        n_chrom = 4
        dd = 0.8 * w
        xs = np.linspace(com[0] - dd * np.cos(chrom_angle), com[0] + dd * np.cos(chrom_angle), n_chrom)
        ys = np.linspace(com[1] - dd * np.sin(chrom_angle), com[1] + dd * np.sin(chrom_angle), n_chrom)

        for i in range(n_chrom):
            line1, line2 = chrom(xs[i], ys[i], chrom_angle)
            plt.plot([spindle_poles[0, 0], xs[i]], [spindle_poles[0, 1], ys[i]], color='g', linewidth=6)
            plt.plot([spindle_poles[1, 0], xs[i]], [spindle_poles[1, 1], ys[i]], color='g', linewidth=6)
            plt.plot(line1[0], line1[1], color='dodgerblue', linewidth=10)
            plt.plot(line2[0], line2[1], color='dodgerblue', linewidth=10)
    else:
        chrom_angle = np.pi / 2 + spindle_angle
        n_chrom = 4
        dd = 0.7 * w
        xs = np.linspace(com[0] - dd * np.cos(chrom_angle), com[0] + dd * np.cos(chrom_angle), n_chrom)
        ys = np.linspace(com[1] - dd * np.sin(chrom_angle), com[1] + dd * np.sin(chrom_angle), n_chrom)

        for i in range(n_chrom):
            line1, line2 = chrom(xs[i], ys[i], chrom_angle)
            plt.plot([spindle_poles[0, 0], xs[i]], [spindle_poles[0, 1], ys[i]], color='g', linewidth=6)
            plt.plot([spindle_poles[1, 0], xs[i]], [spindle_poles[1, 1], ys[i]], color='g', linewidth=6)
            plt.plot(line1[0], line1[1], color='dodgerblue', linewidth=10)
            plt.plot(line2[0], line2[1], color='dodgerblue', linewidth=10)
    ax.legend(ncol=1, loc='upper left', facecolor='white', framealpha=1)

    # SAVING
    path = params[5]
    name = params[4]
    plt.savefig(
        path + '/' + name + '_' + str(1 + n) + '.pdf',
        bbox_inches=None,  # Disable auto-bounding box
        pad_inches=0  # No padding
    )
    plt.close(fig)


# ============================================================
# Main simulation loop — cell-type specific (differs in initial plot_list
# seeding, severance handling, save-interval cadence, and which task_id
# gate is used for PDF saving)
# ============================================================

def simulate_follicle(params):
    # before simulating create cell

    cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, spots, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2 = make_cell_follicle(
        params)
    borders = [min(cell[:, 0]) - 2, max(cell[:, 0]) + 2, max(cell[:, 1]) + 1, min(cell[:, 1]) - 1]
    v_c = np.array([0, 0])
    push_force, pull_force = find_force(astral_MTs, spindle_poles, which_bind, which_push, v_c)
    # plot intial cell and spindle
    if SAVE_FRAME_PLOTS:
        plot_cell_follicle(cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, push_force, pull_force,
                           spots, -1, which_bind, which_push, v_c, -1, 0, params, borders)
    # Dataframe
    df = pd.DataFrame(df_list2,
                      columns=['Run', 'Number', 'Death', 'Switch', 'short', 'bind status', 'Push', 'State', 'Force',
                               'Length', 'Orig length', 'end if out', 'astral angle'])
    # CUTOFF
    plot_list, spindle_center, RMS_angle_list, RMS_center_list = ([] for i in range(4))
    i = 0
    ratio_list = []
    cutoff = 0

    while (1 > 0):
        print(f'{i}-----------------------=++++++++++++++++++++=')
        if (i > int(total_time / time_step)):
            cutoff = i
            break
        # MOVE SPINDLE
        new_cell, new_spots, new_poles, new_angle, new_astral_MTs, new_astral_angles, new_state, new_which_push, new_which_bind, new_free_spots, new_astral_which_spot, new_orig_length, df_list2, ratio, push_force, pull_force, new_v_c = move_spindle(
            params, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, cell, spots, which_push, which_bind,
            free_spots, astral_which_spot, orig_length, v_c, i)

        new_center = np.array(
            [1 / 2 * (new_poles[0, 0] + new_poles[1, 0]), 1 / 2 * (new_poles[0, 1] + new_poles[1, 1])])

        pull_t = np.sum(pull_force[0] + pull_force[1], axis=0)
        push_t = np.sum(push_force[0] + push_force[1], axis=0)

        for pole in new_poles:
            cell_polygon = Polygon(cell)
            point = pole
            if not cell_polygon.contains(Point(point)):
                break
                return plot_list, df

        # DATA LIST
        plot_list.append([cell, new_astral_MTs, new_state, new_poles, new_angle, new_center, new_center[1], spots,
                          np.sum(new_which_bind == 1, axis=1), np.sum(new_which_push == 1, axis=1), 0, ratio, pull_t,
                          push_t])

        # PLOTTING
        save_interval_time = 2  # Save every 0.5 time units
        save_interval_steps = int(round(save_interval_time / time_step))  # Steps between saves
        if i % save_interval_steps == 0 and task_id < 10 and SAVE_FRAME_PLOTS:  # Only save for task_id < 10
            plot_cell_follicle(new_cell, new_astral_MTs, new_astral_angles, new_state, new_poles, new_angle, push_force,
                               pull_force, new_spots, i, new_which_bind, new_which_push, v_c, i, ratio, params, borders)

        cell, spots, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, which_push, which_bind, free_spots, astral_which_spot, orig_length, v_c = new_cell, new_spots, new_astral_MTs, new_astral_angles, new_state, new_poles, new_angle, new_which_push, new_which_bind, new_free_spots, new_astral_which_spot, new_orig_length, new_v_c

        # Dataframes
        temporary_df = pd.DataFrame(df_list2,
                                    columns=['Run', 'Number', 'Death', 'Switch', 'short', 'bind status', 'Push',
                                             'State', 'Force', 'Length', 'Orig length', 'end if out', 'astral angle'])
        df = pd.concat([df, temporary_df])

        i += 1

    return plot_list, df


def simulate_celegans(params):
    # before simulating create cell

    cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, spots, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2 = make_cell_celegans(
        params)
    borders = [min(cell[:, 0]) - 2, max(cell[:, 0]) + 2, max(cell[:, 1]) + 1, min(cell[:, 1]) - 1]
    v_c = np.array([0, 0])
    push_force, pull_force = find_force(astral_MTs, spindle_poles, which_bind, which_push, v_c)
    # plot intial cell and spindle
    if SAVE_FRAME_PLOTS:
        plot_cell_celegans(cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, push_force, pull_force,
                           spots, -1, which_bind, which_push, v_c, -1, 0, params, borders)
    # Dataframe
    df = pd.DataFrame(df_list2,
                      columns=['Run', 'Number', 'Death', 'Switch', 'short', 'bind status', 'Push', 'State', 'Force',
                               'Length', 'Orig length', 'end if out', 'astral angle'])
    # CUTOFF
    plot_list, spindle_center, RMS_angle_list, RMS_center_list = ([] for i in range(4))
    i = 0
    ratio_list = []
    cutoff = 0
    center = np.array(
        [1 / 2 * (spindle_poles[0, 0] + spindle_poles[1, 0]), 1 / 2 * (spindle_poles[0, 1] + spindle_poles[1, 1])])
    plot_list.append([cell, astral_MTs, state, spindle_poles, spindle_angle, center, center[1], spots,
                      np.sum(which_bind == 1, axis=1), np.sum(which_push == 1, axis=1), 0, 0, np.array([0, 0]),
                      np.array([0, 0])])

    while (1 > 0):
        if (i > int(total_time / time_step)):
            cutoff = i
            break
        # MOVE SPINDLE
        if (i * time_step > severance_time and severance == 1):
            new_cell, new_spots, new_poles, new_angle, new_astral_MTs, new_astral_angles, new_state, new_which_push, new_which_bind, new_free_spots, new_astral_which_spot, new_orig_length, df_list2, ratio, push_force, pull_force, new_v_c = move_severed_spindle(
                params, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, cell, spots, which_push,
                which_bind, free_spots, astral_which_spot, orig_length, v_c, i)
        else:
            new_cell, new_spots, new_poles, new_angle, new_astral_MTs, new_astral_angles, new_state, new_which_push, new_which_bind, new_free_spots, new_astral_which_spot, new_orig_length, df_list2, ratio, push_force, pull_force, new_v_c = move_spindle(
                params, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, cell, spots, which_push,
                which_bind, free_spots, astral_which_spot, orig_length, v_c, i)

        new_center = np.array(
            [1 / 2 * (new_poles[0, 0] + new_poles[1, 0]), 1 / 2 * (new_poles[0, 1] + new_poles[1, 1])])

        pull_t = np.sum(pull_force[0] + pull_force[1], axis=0)
        push_t = np.sum(push_force[0] + push_force[1], axis=0)

        for pole in new_poles:
            cell_polygon = Polygon(cell)
            point = pole
            if not cell_polygon.contains(Point(point)):
                break
                return plot_list, df

        # DATA LIST
        if (severance == 1):
            plot_list.append([cell, new_astral_MTs, new_state, new_poles, new_angle, new_poles[0], new_poles[1], spots,
                              np.sum(new_which_bind == 1, axis=1), np.sum(new_which_push == 1, axis=1), 0, ratio,
                              pull_t, push_t])
        else:
            plot_list.append([cell, new_astral_MTs, new_state, new_poles, new_angle, new_center, new_center[1], spots,
                              np.sum(new_which_bind == 1, axis=1), np.sum(new_which_push == 1, axis=1), 0, ratio,
                              pull_t, push_t])

        # PLOTTING
        save_interval_time = 1  # Save every 0.5 time units
        save_interval_steps = int(round(save_interval_time / time_step))  # Steps between saves
        if (i % save_interval_steps == 0 and task_id < 10 and SAVE_FRAME_PLOTS):
            plot_cell_celegans(new_cell, new_astral_MTs, new_astral_angles, new_state, new_poles, new_angle, push_force,
                               pull_force, new_spots, i, new_which_bind, new_which_push, v_c, i, ratio, params, borders)

        cell, spots, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, which_push, which_bind, free_spots, astral_which_spot, orig_length, v_c = new_cell, new_spots, new_astral_MTs, new_astral_angles, new_state, new_poles, new_angle, new_which_push, new_which_bind, new_free_spots, new_astral_which_spot, new_orig_length, new_v_c

        # Dataframes
        temporary_df = pd.DataFrame(df_list2,
                                    columns=['Run', 'Number', 'Death', 'Switch', 'short', 'bind status', 'Push',
                                             'State', 'Force', 'Length', 'Orig length', 'end if out', 'astral angle'])
        df = pd.concat([df, temporary_df])

        i += 1

    return plot_list, df


def simulate_zebrafish(params):
    # before simulating create cell

    cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, spots, which_push, which_bind, free_spots, astral_which_spot, orig_length, df_list2 = make_cell_zebrafish(
        params)
    borders = [min(cell[:, 0]) - 2, max(cell[:, 0]) + 2, max(cell[:, 1]) + 1, min(cell[:, 1]) - 1]
    v_c = np.array([0, 0])
    push_force, pull_force = find_force(astral_MTs, spindle_poles, which_bind, which_push, v_c)
    # plot intial cell and spindle
    if SAVE_FRAME_PLOTS:
        plot_cell_zebrafish(cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, push_force, pull_force,
                            spots, -1, which_bind, which_push, v_c, -1, 0, params, borders)
    # Dataframe
    df = pd.DataFrame(df_list2,
                      columns=['Run', 'Number', 'Death', 'Switch', 'short', 'bind status', 'Push', 'State', 'Force',
                               'Length', 'Orig length', 'end if out', 'astral angle'])
    # CUTOFF
    plot_list, spindle_center, RMS_angle_list, RMS_center_list = ([] for i in range(4))
    i = 0
    ratio_list = []
    cutoff = 0

    while (1 > 0):
        print(f'{i}-----------------------=++++++++++++++++++++=')
        if (i > int(total_time / time_step)):
            cutoff = i
            break
        # MOVE SPINDLE
        new_cell, new_spots, new_poles, new_angle, new_astral_MTs, new_astral_angles, new_state, new_which_push, new_which_bind, new_free_spots, new_astral_which_spot, new_orig_length, df_list2, ratio, push_force, pull_force, new_v_c = move_spindle(
            params, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, cell, spots, which_push, which_bind,
            free_spots, astral_which_spot, orig_length, v_c, i)

        new_center = np.array(
            [1 / 2 * (new_poles[0, 0] + new_poles[1, 0]), 1 / 2 * (new_poles[0, 1] + new_poles[1, 1])])

        pull_t = np.sum(pull_force[0] + pull_force[1], axis=0)
        push_t = np.sum(push_force[0] + push_force[1], axis=0)

        # DATA LIST
        plot_list.append([cell, new_astral_MTs, new_state, new_poles, new_angle, new_center, new_center[1], spots,
                          np.sum(new_which_bind == 1, axis=1), np.sum(new_which_push == 1, axis=1), 0, ratio, pull_t,
                          push_t])

        # PLOTTING
        save_interval_time = 2  # Save every 0.5 time units
        save_interval_steps = int(round(save_interval_time / time_step))  # Steps between saves
        if (i % save_interval_steps == 0 and task_id == 1 and SAVE_FRAME_PLOTS):
            plot_cell_zebrafish(new_cell, new_astral_MTs, new_astral_angles, new_state, new_poles, new_angle,
                                push_force, pull_force, new_spots, i, new_which_bind, new_which_push, v_c, i, ratio,
                                params, borders)

        cell, spots, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, which_push, which_bind, free_spots, astral_which_spot, orig_length, v_c = new_cell, new_spots, new_astral_MTs, new_astral_angles, new_state, new_poles, new_angle, new_which_push, new_which_bind, new_free_spots, new_astral_which_spot, new_orig_length, new_v_c

        # Dataframes
        temporary_df = pd.DataFrame(df_list2,
                                    columns=['Run', 'Number', 'Death', 'Switch', 'short', 'bind status', 'Push',
                                             'State', 'Force', 'Length', 'Orig length', 'end if out', 'astral angle'])
        df = pd.concat([df, temporary_df])

        i += 1

    return plot_list, df


# ============================================================
# Top-level simulate+extract driver — cell-type specific (differs in
# whether/how the raw dataframe is exported and whether angles are flipped)
# ============================================================

def plot_simulate_follicle(params):
    pl_list, df = simulate_follicle(params)
    if task_id == 1 and SAVE_RAW_MT_LOG:
        file_name = params[4] + '.xlsx'
        file_path = os.path.join(params[5], file_name)
        df.head(10000).to_excel(file_path, index=False)

    # EXTRACTING DATA
    df_angle = [np.rad2deg(row[4]) for row in pl_list]
    df_angle = flip_angles(df_angle)
    df_center = [row[5] for row in pl_list]
    df_ypos = [row[6] for row in pl_list]
    df_count_push = [np.sum(row[9]) for row in pl_list]
    df_count_pull = [np.sum(row[8]) for row in pl_list]
    df_ratio = [row[11] for row in pl_list]
    df_pull_t = [row[12] for row in pl_list]
    df_push_t = [row[13] for row in pl_list]
    return df_angle, df_count_push, df_count_pull, df_center, df_ypos, df_ratio, df_pull_t, df_push_t


def plot_simulate_celegans(params):
    pl_list, df = simulate_celegans(params)
    if task_id == 1 and SAVE_RAW_MT_LOG:
        file_name = params[4] + '.xlsx'
        file_path = os.path.join(params[5], file_name)

    # EXTRACTING DATA
    df_angle = [np.rad2deg(row[4]) for row in pl_list]
    df_angle = flip_angles(df_angle)
    df_center = [row[5] for row in pl_list]
    df_ypos = [row[6] for row in pl_list]
    df_count_push = [np.sum(row[9]) for row in pl_list]
    df_count_pull = [np.sum(row[8]) for row in pl_list]
    df_ratio = [row[11] for row in pl_list]
    df_pull_t = [row[12] for row in pl_list]
    df_push_t = [row[13] for row in pl_list]
    return df_angle, df_count_push, df_count_pull, df_center, df_ypos, df_ratio, df_pull_t, df_push_t


def plot_simulate_zebrafish(params):
    pl_list, df = simulate_zebrafish(params)

    # EXTRACTING DATA
    df_angle = [np.rad2deg(row[4]) for row in pl_list]
    df_center = [row[5] for row in pl_list]
    df_ypos = [row[6] for row in pl_list]
    df_count_push = [np.sum(row[9]) for row in pl_list]
    df_count_pull = [np.sum(row[8]) for row in pl_list]
    df_ratio = [row[11] for row in pl_list]
    df_pull_t = [row[12] for row in pl_list]
    df_push_t = [row[13] for row in pl_list]
    return df_angle, df_count_push, df_count_pull, df_center, df_ypos, df_ratio, df_pull_t, df_push_t


# ============================================================
# PUBLIC ENTRY POINT — run a single simulation from explicit parameters
#
# Everything below this point in the __main__ block used to set up module
# globals (cell_type, a, b, r, w, spread, config, the catastrophe/rescue/
# dynamic-binding probabilities, ...) directly from regex tokens pulled out
# of a Slurm-array-friendly run-name string. That convention is kept for the
# CLI (see __main__ below) since existing job-array submissions rely on it,
# but it's not something a new user should have to reverse-engineer just to
# run one simulation.
#
# run_simulation() is the same setup + run logic, factored out so it can be
# called with plain, named parameters instead. It sets the same module
# globals the rest of this file already expects (which is why it must live
# in this module rather than importing this file's functions elsewhere), and
# is what both __main__ and Spindle_Simulator.ipynb call into. Nothing about
# the physics itself changes here — this only changes how a run is configured.
# ============================================================

def run_simulation(
    cell_type,
    *,
    neuroblast=False,             # FE only: True -> apical-only neuroblast FG placement
    celegans_mode="PNC",          # "PNC" | "spindle" (celegans only)
    endo_mode="junc",             # "junc" | "ever" (endo only)
    endo_cell_number=None,        # required for cell_type="endo"; see list_available_endo_cells()
    motor_density=10,
    n_astral_mt=50,
    push=1,
    pull=0,
    spindle_length=None,           # SL, in model units; FE and celegans "spindle" mode only.
                                   # None -> cell-type default (0.8 for FE == 8 µm; 1.8 for celegans "spindle").
                                   # 1 model unit == 10 µm throughout this file.
    initial_angle=None,           # degrees; None -> cell-type default (45 FE, 90 PNC, 0 celegans-spindle;
                                   # ignored for endo, which reads the tracked initial angle from data/Movie_info.xlsx)
    time_step=0.05,
    RC_rate=1.0,                  # catastrophe/rescue rate multiplier (FE/endo only; celegans hardcodes rates)
    GS_rate=1.0,                  # growth/shrink rate multiplier (FE/endo only)
    total_time=None,              # None -> cell-type default (matches the original per-branch literals)
    advanced=None,                # dict of SimulationParameters attribute overrides (mu_fric, visc, EI, ...)
    seed=None,
    run_label="run",
    test_folder_path="test_out",
    data_folder_path="data_out",
    task_id=1,
    save_frame_plots=True,
    save_raw_mt_log=True,
    verbose=True,
):
    """
    Run one simulation and return its result dataframes as a dict.

    This is the explicit-parameter equivalent of the CLI's run-name string —
    see Spindle_Simulator.ipynb for a widget-driven interface built on top of
    this function, and CLAUDE.md for the parameter reference.
    """
    if cell_type not in ("FE", "celegans", "endo"):
        raise ValueError(f"Unknown cell_type: {cell_type!r} (expected 'FE', 'celegans', or 'endo')")
    if cell_type == "endo" and endo_cell_number is None:
        raise ValueError("endo_cell_number is required for cell_type='endo' (see list_available_endo_cells())")

    g = globals()
    g['SAVE_FRAME_PLOTS'] = save_frame_plots
    g['SAVE_RAW_MT_LOG'] = save_raw_mt_log
    g['task_id'] = task_id
    g['cell_type'] = cell_type
    g['push'] = push
    g['pull'] = pull
    g['time_step'] = time_step
    # Always defined regardless of cell type: SimulationParameters.__init__
    # reads both RC_rate and GS_rate together for its 'FE'/'endo' branch, but
    # the original run-name parsers only ever return one or the other
    # (extract_parameters_follicle -> RC_rate only, extract_parameters_zebrafish
    # -> GS_rate only) — leaving the other undefined would raise a NameError.
    # Defaulting both to 1.0 here keeps every path safe.
    g['RC_rate'] = RC_rate
    g['GS_rate'] = GS_rate

    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    job_id = os.getenv('SLURM_JOB_ID')

    config = SimulationParameters()
    g['config'] = config
    if advanced:
        for key, value in advanced.items():
            setattr(config, key, value)

    prob_catastr = 1 - np.exp(-config.catastr_rate * time_step)
    prob_rescue = 1 - np.exp(-config.rescue_rate * time_step)
    if pull != 0:
        prob_dyn_bind = 1 - np.exp(-config.dyn_bind * time_step)
        prob_dyn_unbind = 1 - np.exp(-config.dyn_unbind * time_step)
    else:
        prob_dyn_bind = 0
        prob_dyn_unbind = 1
    g['prob_catastr'] = prob_catastr
    g['prob_rescue'] = prob_rescue
    g['prob_dyn_bind'] = prob_dyn_bind
    g['prob_dyn_unbind'] = prob_dyn_unbind

    if cell_type == "celegans":
        model = celegans_mode
        g['model'] = model
        stall_time = 0
        g['stall_time'] = stall_time
        g['severance_time'] = 12222220
        g['severance'] = 0
        g['mu'] = config.mu_fric

        if model == 'PNC':
            default_angle = 90
            g['total_time'] = total_time if total_time is not None else 360 + stall_time
            a, b, r, w = 2.5, 1.5, 0.5, 0.5
            g['spread'] = np.deg2rad(180)
        elif model == 'spindle':
            default_angle = 0
            SL = spindle_length if spindle_length is not None else 1.8
            r = SL / 2
            g['total_time'] = total_time if total_time is not None else 300 + stall_time
            a, b, w = 2.5, 1.5, 0.5
            g['spread'] = np.deg2rad(270)
        else:
            raise ValueError(f"Unknown celegans_mode: {celegans_mode!r} (expected 'PNC' or 'spindle')")
        g['a'], g['b'], g['r'], g['w'] = a, b, r, w
        g['mean'], g['stdev'] = 1, 0.3
        g['AL'], g['scale'] = 9, 0.167
        angle = initial_angle if initial_angle is not None else default_angle
        sp_angle = np.deg2rad(angle) * np.ones(50)
        extranote = "make astral check bind init"

    elif cell_type == "endo":
        c = endo_cell_number
        g['model'] = endo_mode
        g['c'] = c
        df_exp = pd.read_excel(os.path.join(DATA_DIR, 'cells_data.xlsx'), sheet_name='cell_' + str(c))
        g['real_exp_init_angle'] = df_exp['Spindle angle'].values[0]
        df_data = pd.read_excel(os.path.join(DATA_DIR, 'Movie_info.xlsx'), sheet_name='Alikhan', index_col=False)
        g['df_data'] = df_data
        frame_rates = df_data['Frame rate']
        ends = df_data['anaphase onset']
        starts = df_data['metaphase']
        rescale = df_data['model scale factor']
        initial_angles = df_data['Initial angle']
        spindle_length_df = df_data['Spindle length']
        spindle_width_df = df_data['Spindle width']
        g['frame_rates'], g['ends'], g['starts'], g['rescale'] = frame_rates, ends, starts, rescale
        g['initial_angles'] = initial_angles
        g['final_angle'] = df_data['Final angle']
        g['spindle_length_df'], g['spindle_width_df'] = spindle_length_df, spindle_width_df
        g['total_time'] = total_time if total_time is not None else (ends[c - 1] - starts[c - 1]) * frame_rates[c - 1]
        g['long_axis'] = df_data['long axis PCA'][c - 1]
        a, b = 1, 1
        g['spread'] = 3 * np.pi / 2
        r = 0.1 * spindle_length_df[c - 1] / 2
        w = 0.1 * spindle_width_df[c - 1] / 2
        g['a'], g['b'], g['r'], g['w'] = a, b, r, w
        g['AL'], g['scale'] = 1, 1
        angle = initial_angle if initial_angle is not None else initial_angles[c - 1]
        sp_angle = np.deg2rad(angle) * np.ones(50)
        extranote = "ever motors new check spindle"

    else:  # "FE"
        g['NB'] = 1 if neuroblast else 0
        g['severance_time'] = 5000
        g['severance'] = 0
        g['mu'] = config.mu_fric
        g['total_time'] = total_time if total_time is not None else 600
        a = b = 0.5
        SL = spindle_length if spindle_length is not None else 0.8
        if SL > 0.95:
            raise ValueError(
                f"Spindle length {SL:.2f} model units ({SL * 10:.1f} µm) is too long for a follicle "
                f"epithelial cell (cell semi-axes a=b={a} model units) — use spindle_length <= 0.95 model "
                f"units (<= 9.5 µm)."
            )
        r = SL / 2
        g['spread'] = 3 * np.pi / 2
        w = 0.8 * a
        g['a'], g['b'], g['r'], g['w'] = a, b, r, w
        g['AL'], g['scale'] = 1, 1
        angle = initial_angle if initial_angle is not None else 90
        sp_angle = np.deg2rad(angle) * np.ones(50)
        extranote = "mean MT length 1"

    folder_name, new_dir_path = create_simulation_directory(test_folder_path, run_label, task_id)

    additional_params = {
        "job_id": job_id,
        "total_time": g['total_time'],
        "a": a, "b": b, "r": r, "spread": g['spread'], "w": w,
        "mean MT length": g['AL'], "stdev MT length": g['scale'],
        "sp_angle": sp_angle[0],
        "push": push,
        "length_MTs": config.length_MTs,
        "state_MTs": config.state_MTs,
        "extranote": extranote,
    }

    def _run_and_save():
        config.save_parameters_to_file(data_folder_path, run_label + "_parameters.txt", additional_params)
        run_params = [[a, b, r, w], motor_density, n_astral_mt, sp_angle[int(task_id - 1)], folder_name, new_dir_path]
        if cell_type == "celegans":
            result = plot_simulate_celegans(run_params)
        elif cell_type == "endo":
            result = plot_simulate_zebrafish(run_params)
        else:
            result = plot_simulate_follicle(run_params)

        df_angle, df_count_push, df_count_pull, df_center, df_ypos, df_ratio, df_pull_t, df_push_t = result
        data_dict = {
            "Angle": df_angle, "Center": df_center, "Y-pos": df_ypos,
            "N_pull": df_count_pull, "N_push": df_count_push, "Ratio": df_ratio,
            "Pull_t": df_pull_t, "Push_t": df_push_t,
        }
        excel_file_path = os.path.join(data_folder_path, f"{folder_name}_data.xlsx")
        with pd.ExcelWriter(excel_file_path, engine="openpyxl") as writer:
            for sheet_name, data in data_dict.items():
                pd.DataFrame({"Run " + str(task_id): data}).to_excel(writer, sheet_name=sheet_name, index=False)
        return data_dict, excel_file_path

    if verbose:
        data_dict, excel_file_path = _run_and_save()
    else:
        with contextlib.redirect_stdout(io.StringIO()):
            data_dict, excel_file_path = _run_and_save()

    return {
        "angle": data_dict["Angle"],
        "count_push": data_dict["N_push"],
        "count_pull": data_dict["N_pull"],
        "center": data_dict["Center"],
        "ypos": data_dict["Y-pos"],
        "ratio": data_dict["Ratio"],
        "pull_t": data_dict["Pull_t"],
        "push_t": data_dict["Push_t"],
        "time_step": time_step,
        "folder_name": folder_name,
        "run_dir": new_dir_path,
        "data_xlsx": excel_file_path,
    }


# ============================================================
# CLI entry point — dispatches to the correct cell-type family based on the
# 'cell_' token embedded in the run name (each original script only ever
# saw its own family's naming convention and run-name parser; this probe
# reproduces that dispatch inside a single unified script), then hands off
# to run_simulation() above.
# ============================================================

if __name__ == "__main__":

    test_folder_path = sys.argv[1]
    data_folder_path = sys.argv[2]
    name = sys.argv[3]  # "spindle_9214_FE_AL_1_0.3_SL_1.8_ts_0.05_opt_forces_3.6_4_2.5_MUD_PINS_10_MT_50_push_1"
    # "spindle_9214_endo_cell_1_AL_1_0.3_SL_1.8_ts_0.05_opt_forces_3.6_4_2.5_MUD_PINS_10_MT_50_push_1"

    print(f"Running simulation with name: {name}")
    # Falls back to task_id=1 / job_id=None when not launched under sbatch, so
    # this script also runs as a plain local command for a single simulation.
    task_id = int(os.getenv('SLURM_ARRAY_TASK_ID', 1))
    job_id = os.getenv('SLURM_JOB_ID')

    # Probe which family's run-name convention applies (this mirrors which of
    # the three original scripts would have been invoked for this run)
    _cell_probe_match = re.search(r'cell_([A-Za-z]+)', name)
    _cell_type_probe = _cell_probe_match.group(1) if _cell_probe_match else None

    if _cell_type_probe == "celegans":
        motor_density, astral_number, push, pull, SL, time_step, cell_type, init_sp_angle, model = extract_parameters_celegans(
            name)
        run_kwargs = dict(cell_type=cell_type, celegans_mode=model, spindle_length=SL, initial_angle=init_sp_angle)
    elif _cell_type_probe == "endo":
        # NB: the parsed angle (if any) is intentionally NOT forwarded here —
        # the original script always used the tracked initial angle from
        # Movie_info.xlsx for endo runs regardless of any angle_ token in the
        # run name, and run_simulation() replicates that by only falling back
        # to the tracked angle when initial_angle is left as None.
        motor_density, astral_number, push, pull, _, time_step, cell_type, init_sp_angle, model, GS_rate, AL_length = extract_parameters_zebrafish(
            name)
        run_kwargs = dict(cell_type=cell_type, endo_mode=model, endo_cell_number=find_number_after_cell(name),
                          GS_rate=GS_rate)
    elif _cell_type_probe == "FE":
        motor_density, astral_number, push, pull, SL, time_step, cell_type, init_sp_angle, NB, RC_rate = extract_parameters_follicle(
            name)
        run_kwargs = dict(cell_type=cell_type, neuroblast=bool(NB), spindle_length=SL, RC_rate=RC_rate,
                          initial_angle=init_sp_angle)
    else:
        raise ValueError(f"Unknown cell type in run name: {name}")

    run_simulation(
        motor_density=motor_density,
        n_astral_mt=astral_number,
        push=push,
        pull=pull,
        time_step=time_step,
        run_label=name,
        test_folder_path=test_folder_path,
        data_folder_path=data_folder_path,
        task_id=task_id,
        # Slurm/batch runs keep producing every artifact, exactly like before
        # this file grew a notebook-friendly API.
        save_frame_plots=True,
        save_raw_mt_log=True,
        verbose=True,
        **run_kwargs,
    )


# ============================================================
# NOTEBOOK / INTERACTIVE UI HELPERS
#
# Everything below is additive and only used by Spindle_Simulator.ipynb (or
# any other interactive front-end someone wants to build on run_simulation()
# above) — none of it is imported or exercised by the CLI/Slurm path. Only
# build_ui() needs ipywidgets, imported lazily so importing this module for
# a batch run never requires it to be installed.
# ============================================================

def list_available_endo_cells():
    """
    Tracked-cell numbers that can actually be run as cell_type="endo": present
    as data/junctions/C<n>_Mask_Movie_Junctions.txt (ignoring macOS '._'-prefixed
    artifact files already sitting in that folder) AND having a 'cell_<n>'
    sheet in data/cells_data.xlsx. Used to populate the notebook's tracked-cell
    dropdown with only choices that will actually run.
    """
    junction_dir = os.path.join(DATA_DIR, "junctions")
    if not os.path.isdir(junction_dir):
        return []
    found = set()
    for fname in os.listdir(junction_dir):
        if fname.startswith("._"):
            continue
        m = re.match(r"C(\d+)_Mask_Movie_Junctions", fname)
        if m:
            found.add(int(m.group(1)))

    cells_data_path = os.path.join(DATA_DIR, "cells_data.xlsx")
    if not found or not os.path.exists(cells_data_path):
        return []
    wb = openpyxl.load_workbook(cells_data_path, read_only=True)
    try:
        available_sheets = set(wb.sheetnames)
    finally:
        wb.close()

    return sorted(c for c in found if f"cell_{c}" in available_sheets)


def plot_summary(results, title=None):
    """
    Plot the summary curves from a run_simulation(...) result dict: spindle
    angle, push/pull microtubule counts, spindle-center x/y position, and the
    2D center trajectory, all vs. time (except the trajectory).
    """
    dt = results["time_step"]
    center = np.asarray(results["center"])
    t = np.arange(len(results["angle"])) * dt

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    if title:
        fig.suptitle(title)

    ax = axes[0, 0]
    ax.plot(t, results["angle"])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Spindle angle (deg)")
    ax.set_title("Spindle angle vs. time")

    ax = axes[0, 1]
    ax.plot(t, results["count_push"], label="Pushing MTs")
    ax.plot(t, results["count_pull"], label="Pulling MTs")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Count")
    ax.set_title("Push / pull microtubule counts vs. time")
    ax.legend()

    ax = axes[1, 0]
    ax.plot(t, center[:, 0], label="x position")
    ax.plot(t, results["ypos"], label="y position")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Position")
    ax.set_title("Spindle center position vs. time")
    ax.legend()

    ax = axes[1, 1]
    ax.plot(center[:, 0], results["ypos"], '-', color="gray", linewidth=1)
    ax.plot(center[0, 0], results["ypos"][0], 'go', label="start")
    ax.plot(center[-1, 0], results["ypos"][-1], 'ro', label="end")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Spindle center trajectory")
    ax.set_aspect('equal', adjustable='datalim')
    ax.legend()

    fig.tight_layout()
    # Deliberately no plt.show(): under a real Jupyter inline backend that's
    # harmless, but outside one (e.g. a plain script, or an unusual backend)
    # it can select an interactive GUI backend and block. Returning the
    # figure lets Jupyter's normal rich-display auto-show it, and lets
    # build_ui() explicitly display()+close() it inside the Output widget.
    return fig


def build_ui():
    """
    Build the widget-driven interface used by Spindle_Simulator.ipynb. Usage
    in a notebook cell:

        ui = build_ui()
        display(ui)
    """
    import ipywidgets as w
    from IPython.display import display

    label_w = {"description_width": "150px"}
    full_w = w.Layout(width="420px")

    CELL_TYPES = {
        "Fly follicular epithelium": ("FE", False),
        "Fly neuroblast": ("FE", True),
        "C. elegans — pronuclei centering": ("celegans", "PNC"),
        "C. elegans — spindle positioning": ("celegans", "spindle"),
        "Zebrafish embryo (tracked cell)": ("endo", None),
    }

    # Rough per-astral-MT-and-cell-type simulation cost, fit from measured
    # wall-clock timings on one machine (ms per integration step = BASE +
    # SLOPE * n_astral_mt). Only meant to give the user a ballpark before they
    # commit to a long run -- actual speed varies with hardware.
    TIMING_FIT = {"FE": (8, 6.4), "celegans": (2.5, 2.0), "endo": (4, 3.2)}

    cell_type_dd = w.Dropdown(options=list(CELL_TYPES), description="Cell type:", style=label_w, layout=full_w)
    angle_slider = w.FloatSlider(value=90, min=0, max=180, step=1, description="Initial angle (°):",
                                 style=label_w, layout=full_w)
    # Follicle epithelium/neuroblast take the spindle length in µm (1 model
    # unit == 10 µm; default 8 µm == 0.8 model units); celegans "spindle" mode
    # takes it directly in model units (default 1.8). _refresh_fields() swaps
    # the description/default below; _on_run_clicked converts FE's µm value
    # back to model units before calling run_simulation().
    sl_box = w.BoundedFloatText(value=8, min=0.1, max=100, step=0.5, description="Spindle length (µm):",
                                style=label_w, layout=full_w)
    endo_cell_dd = w.Dropdown(options=list_available_endo_cells(), description="Tracked cell #:",
                              style=label_w, layout=full_w)
    endo_info_label = w.HTML(value="")

    # "Number of motors" (FE/celegans) vs "Motor density" (endo) is the same
    # underlying parameter/fill-in box, just interpreted differently downstream:
    # make_fgs_follicle/_neuroblast/_celegans_pnc/_spindle treat it as a plain
    # discrete motor count, while make_fgs_zebrafish_junc/_uf multiply it by
    # the (real, tracked) cell perimeter to get an actual areal/line density.
    # The box's description is switched per cell type below.
    motor_box = w.BoundedFloatText(value=10, min=0, max=100000, step=1, description="Number of motors:",
                                   style=label_w, layout=full_w)
    mt_box = w.BoundedIntText(value=50, min=1, max=100000, step=1, description="Astral microtubules:",
                              style=label_w, layout=full_w)
    duration_box = w.BoundedFloatText(value=600, min=1, max=100000, step=10, description="Duration (s):",
                                      style=label_w, layout=full_w)
    runtime_label = w.HTML(value="")
    push_cb = w.Checkbox(value=True, description="Cortical pushing (astral MTs)")
    pull_cb = w.Checkbox(value=False, description="Cortical pulling (motors)")
    save_outputs_cb = w.Checkbox(value=False,
                                 description="Also save per-frame PDF plots + raw MT log to disk (slower)")
    seed_use_cb = w.Checkbox(value=False, description="Fix random seed")
    seed_box = w.IntText(value=0, description="Seed:", style=label_w, layout=full_w)
    label_box = w.Text(value="my_run", description="Run label:", style=label_w, layout=full_w)
    outdir_box = w.Text(value="notebook_runs", description="Output folder:", style=label_w, layout=full_w)

    rc_slider = w.FloatSlider(value=1.0, min=0.2, max=3.0, step=0.1, description="Catastrophe/rescue rate ×:",
                              style=label_w, layout=full_w)
    gs_slider = w.FloatSlider(value=1.0, min=0.2, max=3.0, step=0.1, description="Growth/shrink rate ×:",
                              style=label_w, layout=full_w)
    mu_fric_box = w.FloatText(value=500, description="Friction (mu_fric):", style=label_w, layout=full_w)
    visc_box = w.FloatText(value=100, description="Viscosity:", style=label_w, layout=full_w)
    ei_box = w.FloatText(value=0.2, description="Rigidity (EI):", style=label_w, layout=full_w)
    advanced = w.Accordion(children=[w.VBox([rc_slider, gs_slider, mu_fric_box, visc_box, ei_box])])
    advanced.set_title(0, "Advanced (physical parameters)")
    advanced.selected_index = None

    basic_box = w.VBox([])
    run_btn = w.Button(description="Run simulation", button_style="success", icon="play")
    status = w.HTML(value="")
    output = w.Output()

    def _refresh_runtime_estimate(*_):
        cell_type, _sub = CELL_TYPES[cell_type_dd.value]
        base, slope = TIMING_FIT[cell_type]
        ms_per_step = base + slope * mt_box.value
        est_s = (duration_box.value / 0.05) * ms_per_step / 1000
        text = f"~{est_s:.0f} s" if est_s < 90 else f"~{est_s / 60:.1f} min"
        runtime_label.value = f"<i>Very rough runtime estimate: {text} (varies by machine). Lower the duration or astral MT count for a faster look.</i>"

    def _refresh_endo_info(*_):
        if endo_cell_dd.value is None:
            endo_info_label.value = "<i>No tracked cells found (need data/ + data/junctions/).</i>"
            return
        try:
            movie_info = pd.read_excel(os.path.join(DATA_DIR, "Movie_info.xlsx"), sheet_name="Alikhan",
                                       index_col=False)
            c = endo_cell_dd.value
            angle = movie_info["Initial angle"][c - 1]
            duration = (movie_info["anaphase onset"][c - 1] - movie_info["metaphase"][c - 1]) * movie_info["Frame rate"][c - 1]
            duration_box.value = float(duration)
            endo_info_label.value = (f"<i>Initial angle: {angle:.1f}° and duration: {duration:.0f}s "
                                     f"(both read from tracked movie data; duration above is still editable)</i>")
        except Exception:
            endo_info_label.value = "<i>Initial angle and duration are read from tracked movie data.</i>"
        _refresh_runtime_estimate()

    def _refresh_fields(*_):
        cell_type, sub = CELL_TYPES[cell_type_dd.value]
        rows = []
        motor_box.description = "Motor density:" if cell_type == "endo" else "Number of motors:"
        if cell_type == "FE":
            sl_box.description = "Spindle length (µm):"
            sl_box.min, sl_box.max, sl_box.step = 0.1, 100, 0.5
            sl_box.value = 8
            rows = [angle_slider, sl_box]
            duration_box.value = 600
        elif cell_type == "celegans":
            angle_slider.value = 90 if sub == "PNC" else 0
            if sub == "spindle":
                sl_box.description = "Spindle length (model units):"
                sl_box.min, sl_box.max, sl_box.step = 0.4, 2.4, 0.1
                sl_box.value = 1.8
            rows = [angle_slider] + ([sl_box] if sub == "spindle" else [])
            duration_box.value = 360 if sub == "PNC" else 300
        elif cell_type == "endo":
            rows = [endo_cell_dd, endo_info_label]
            _refresh_endo_info()
        rc_slider.disabled = gs_slider.disabled = (cell_type == "celegans")
        basic_box.children = rows
        _refresh_runtime_estimate()

    cell_type_dd.observe(_refresh_fields, names="value")
    endo_cell_dd.observe(_refresh_endo_info, names="value")
    mt_box.observe(_refresh_runtime_estimate, names="value")
    duration_box.observe(_refresh_runtime_estimate, names="value")
    _refresh_fields()

    def _on_run_clicked(_):
        run_btn.disabled = True
        status.value = "<b>Running…</b> see the runtime estimate above; long durations/MT counts can take a while."
        output.clear_output()
        try:
            cell_type, sub = CELL_TYPES[cell_type_dd.value]
            test_folder_path = os.path.join(outdir_box.value, "frames")
            data_folder_path = os.path.join(outdir_box.value, "data")
            os.makedirs(test_folder_path, exist_ok=True)
            os.makedirs(data_folder_path, exist_ok=True)

            kwargs = dict(
                cell_type=cell_type,
                motor_density=motor_box.value,
                n_astral_mt=mt_box.value,
                push=1 if push_cb.value else 0,
                pull=1 if pull_cb.value else 0,
                total_time=duration_box.value,
                run_label=label_box.value or "run",
                test_folder_path=test_folder_path,
                data_folder_path=data_folder_path,
                save_frame_plots=save_outputs_cb.value,
                save_raw_mt_log=save_outputs_cb.value,
                verbose=False,
                seed=(seed_box.value if seed_use_cb.value else None),
                advanced={"mu_fric": mu_fric_box.value, "visc": visc_box.value, "EI": ei_box.value},
            )
            if cell_type == "FE":
                # sl_box is in µm for FE; run_simulation expects model units (1 model unit == 10 µm).
                kwargs.update(neuroblast=sub, spindle_length=sl_box.value / 10.0, initial_angle=angle_slider.value,
                              RC_rate=rc_slider.value, GS_rate=gs_slider.value)
            elif cell_type == "celegans":
                # sl_box is already in model units here.
                kwargs.update(celegans_mode=sub, initial_angle=angle_slider.value,
                              spindle_length=sl_box.value if sub == "spindle" else 1.8)
            else:  # endo
                if endo_cell_dd.value is None:
                    raise ValueError("No tracked cell available — check that data/ and data/junctions/ are present.")
                kwargs.update(endo_cell_number=endo_cell_dd.value, GS_rate=gs_slider.value)

            results = run_simulation(**kwargs)
            with output:
                fig = plot_summary(results, title=cell_type_dd.value)
                display(fig)
                plt.close(fig)
            status.value = (f"<b>Done.</b> Data file: <code>{results['data_xlsx']}</code>"
                            + (f" — frame plots in <code>{results['run_dir']}</code>" if save_outputs_cb.value else ""))
        except Exception as exc:
            status.value = f"<span style='color:#b00020'><b>Run failed:</b> {exc}</span>"
        finally:
            run_btn.disabled = False

    run_btn.on_click(_on_run_clicked)

    return w.VBox([
        w.HTML("<h3>1. Choose a cell type</h3>"),
        cell_type_dd, basic_box,
        w.HTML("<h3>2. Simulation settings</h3>"),
        motor_box, mt_box, duration_box, runtime_label, push_cb, pull_cb, save_outputs_cb,
        w.HBox([seed_use_cb, seed_box]),
        label_box, outdir_box,
        advanced,
        run_btn, status, output,
    ])
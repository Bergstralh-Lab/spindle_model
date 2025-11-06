# Spindle Positioning Simulation - User Guide

A comprehensive guide for using the spindle positioning simulation model.

## Table of Contents

1. [Installation](#installation)
2. [Basic Usage](#basic-usage)
3. [Configuration Guide](#configuration-guide)
4. [Cell Types](#cell-types)
5. [Parameters Reference](#parameters-reference)
6. [Running Simulations](#running-simulations)
7. [Output and Analysis](#output-and-analysis)
8. [Advanced Usage](#advanced-usage)
9. [Troubleshooting](#troubleshooting)

---

## Installation

### Requirements

- Python 3.9 or higher
- Required packages:
  - numpy
  - pandas
  - openpyxl
  - matplotlib
  - opencv-python (cv2)
  - shapely
  - scipy

### Setup

```bash
# Install required packages
pip install numpy pandas openpyxl matplotlib opencv-python shapely scipy

# Download the script
# Place spindle.py in your working directory
```

### Verify Installation

```python
python -c "from spindle import simulate; print('Installation successful!')"
```

---

## Basic Usage

### Quickest Start - Run an Example

```python
from spindle import example_celegans_spindle

# Run simulation with default parameters
results = example_celegans_spindle()
```

This will:
- Create a simulation with default parameters
- Run for 10 seconds (simulated time)
- Save plots to `./output/SM_celegans_spindle_run_1/`
- Save data to Excel file

### Simple Custom Simulation

```python
from spindle import CElegansSpindleConfig, simulate

# Create configuration
config = CElegansSpindleConfig(
    time_step=0.05,        # 50 ms timesteps
    total_time=10.0,       # 10 second simulation
    n_astral_mts=100,      # 100 MTs per pole
    pull_force=5.0,        # 5 pN pulling force
    push_force=5.0         # 5 pN pushing force
)

# Run simulation
results = simulate(config, run_id=1)

# Access results
print(f"Final angle: {results['angles'][-1]:.1f}°")
print(f"Final position: {results['centers'][-1]}")
```

---

## Configuration Guide

### Available Cell Types

The model supports five cell type configurations:

1. **FollicleEpithelialConfig** - *Drosophila* follicular epithelium
2. **NeuroblastConfig** - *Drosophila* neuroblast
3. **CElegansPNCConfig** - *C. elegans* pronuclear complex
4. **CElegansSpindleConfig** - *C. elegans* spindle positioning
5. **ZebrafishEndoConfig** - Zebrafish endothelial cell

### Basic Configuration Structure

All configurations share common parameters:

```python
config = SomeCellConfig(
    # Time parameters
    time_step=0.05,           # Simulation timestep (seconds)
    total_time=10.0,          # Total simulation time (seconds)
    seed=42,                  # Random seed for reproducibility
    
    # Spindle geometry
    spindle_half_length=0.4,  # Half-length of spindle (μm)
    spindle_width=0.32,       # Width of spindle (μm)
    spindle_angle_init=0.0,   # Initial angle (radians)
    
    # Microtubules
    n_astral_mts=100,         # Number of MTs per pole
    mt_mean_length=1.0,       # Mean MT length (μm)
    mt_length_stdev=1.0,      # Standard deviation of length
    
    # Forces
    pull_force=5.0,           # Pulling force magnitude (pN)
    push_force=0.0,           # Pushing force magnitude (pN)
    
    # Rates (per second)
    growth_rate=0.013,        # MT growth rate (μm/s)
    shrink_rate=-0.027,       # MT shrinking rate (μm/s)
    catastrophe_rate=0.021,   # Catastrophe frequency (1/s)
    rescue_rate=0.029,        # Rescue frequency (1/s)
    bind_rate=0.03,           # Motor binding rate (1/s)
    unbind_rate=0.02,         # Motor unbinding rate (1/s)
    
    # Physical constants
    viscosity=100.0,          # Cytoplasmic viscosity (Pa·s)
    bending_rigidity=2.0,     # MT bending rigidity (pN·μm²)
    
    # Output
    output_dir=Path("./output"),
    save_plots=True,
    plot_interval=1           # Save every N steps
)
```

---

## Cell Types

### 1. Follicular Epithelial Cell

Models spindle positioning in *Drosophila* egg chamber follicle cells with basal/lateral force generators.

```python
from spindle import FollicleEpithelialConfig, simulate

config = FollicleEpithelialConfig(
    time_step=0.05,
    total_time=30.0,
    
    # Cell geometry (circular)
    cell_radius_a=0.5,        # Cell radius X (μm)
    cell_radius_b=0.5,        # Cell radius Y (μm)
    
    # Force generators
    fg_density=100,           # Number of cortical binding sites
    
    # Spindle
    spindle_half_length=0.4,
    spindle_width=0.32,
    spindle_angle_init=1.57,  # 90° (vertical)
    
    # Microtubules
    n_astral_mts=100,
    mt_mean_length=1.0,
    mt_length_stdev=1.0,
    
    # Forces
    pull_force=5.0,
    push_force=0.0,
    
    output_dir=Path("./output/follicle")
)

results = simulate(config, run_id=1)
```

**Key Features:**
- Gaussian-weighted FG distribution on basal/lateral cortex
- Circular cell geometry
- FGs concentrated around cell equator (±60° from horizontal)

### 2. Neuroblast

Models *Drosophila* neuroblast with apical force generator enrichment.

```python
from spindle import NeuroblastConfig, simulate

config = NeuroblastConfig(
    time_step=0.05,
    total_time=30.0,
    
    # Cell geometry
    cell_radius_a=0.5,
    cell_radius_b=0.5,
    
    # Force generators (separate apical and basal)
    fg_density_basal=100,     # Basal/lateral FGs
    fg_density_apical=100,    # Apical FGs
    
    # Spindle
    spindle_half_length=0.4,
    spindle_width=0.32,
    spindle_angle_init=1.57,  # Start vertical
    
    # Microtubules
    n_astral_mts=100,
    
    # Forces
    pull_force=5.0,
    push_force=0.0,
    
    output_dir=Path("./output/neuroblast")
)

results = simulate(config, run_id=1)
```

**Key Features:**
- Apical FG enrichment (45°-135° from horizontal)
- Models asymmetric division
- Can test apical vs basal force balance

### 3. C. elegans Pronuclear Complex (PNC)

Models early *C. elegans* embryo with posterior-enriched force generators during pronuclear migration.

```python
from spindle import CElegansPNCConfig, simulate

config = CElegansPNCConfig(
    time_step=0.05,
    total_time=10.0,
    
    # Cell geometry (superellipse)
    cell_radius_a=2.5,        # Major axis (μm)
    cell_radius_b=1.5,        # Minor axis (μm)
    superellipse_n=2.2,       # Shape parameter
    
    # Force generators (asymmetric)
    fg_density_anterior=40,   # Fewer FGs anteriorly
    fg_density_posterior=60,  # More FGs posteriorly
    
    # Spindle
    spindle_half_length=0.5,
    spindle_width=0.5,
    # spindle_angle_init auto-set to π/2
    
    # Microtubules (longer in C. elegans)
    n_astral_mts=100,
    mt_mean_length=9.0,
    mt_length_stdev=0.167,
    
    # Forces
    pull_force=5.0,
    push_force=0.0,
    
    # Rates (different from other cell types)
    catastrophe_rate=0.014,
    rescue_rate=0.044,
    
    output_dir=Path("./output/celegans_pnc")
)

results = simulate(config, run_id=1)
```

**Key Features:**
- Elongated superellipse geometry
- Posterior enrichment of force generators
- Longer microtubules
- Different dynamics rates
- Models PNC migration phase

### 4. C. elegans Spindle Positioning

Models *C. elegans* one-cell embryo spindle positioning with uniform FG distribution.

```python
from spindle import CElegansSpindleConfig, simulate

config = CElegansSpindleConfig(
    time_step=0.05,
    total_time=10.0,
    
    # Cell geometry
    cell_radius_a=2.5,
    cell_radius_b=1.5,
    superellipse_n=2.2,
    
    # Force generators (more uniform)
    fg_density_anterior=40,
    fg_density_posterior=60,
    
    # Spindle
    spindle_half_length=0.9,  # Longer spindle
    spindle_width=0.5,
    
    # Microtubules
    n_astral_mts=100,
    mt_mean_length=9.0,
    mt_length_stdev=0.167,
    
    # Forces (both pull and push)
    pull_force=5.0,
    push_force=5.0,           # Pushing active
    
    # Rates
    catastrophe_rate=0.014,
    rescue_rate=0.044,
    
    output_dir=Path("./output/celegans_spindle")
)

results = simulate(config, run_id=1)
```

**Key Features:**
- Similar to PNC but with longer spindle
- Both pulling and pushing forces active
- More uniform FG distribution
- Models metaphase/anaphase positioning

### 5. Zebrafish Endothelial Cell

Models spindle positioning using live imaging data from zebrafish endothelial cells.

```python
from spindle import ZebrafishEndoConfig, simulate
from pathlib import Path

config = ZebrafishEndoConfig(
    time_step=0.05,
    total_time=100.0,         # Longer simulation
    
    # Cell selection
    endo_index=1,             # Which cell (1-15 available)
    
    # Force generator distribution
    fg_distribution="uniform", # "uniform" or "junctions"
    fg_density=10,            # FGs per μm of perimeter
    
    # Microtubules
    n_astral_mts=100,
    mt_mean_length=1.0,
    mt_length_stdev=1.0,
    
    # Forces
    pull_force=5.0,
    push_force=0.0,
    
    # Data paths (adjust to your setup)
    movie_info_file=Path("./data/Movie_info.xlsx"),
    cell_images_dir=Path("./data/cells"),
    spindle_images_dir=Path("./data/spindles"),
    junction_data_dir=Path("./data/junctions"),  # Only needed for junctions mode
    
    output_dir=Path("./output/zebrafish")
)

results = simulate(config, run_id=1)
```

**Key Features:**
- Loads real cell shapes from microscopy images
- Cell shape updates dynamically during simulation
- Spindle geometry from experimental measurements
- Two FG modes: uniform or junction-localized
- Requires experimental data files

**Data Structure Required:**
```
data/
├── Movie_info.xlsx              # Metadata for all cells
├── cells/
│   └── cell_1/
│       ├── Mask_1.jpg
│       ├── Mask_2.jpg
│       └── ...
├── spindles/
│   └── spindle_1/
│       ├── Spindle_1.jpg
│       └── ...
└── junctions/                   # Optional, for junction mode
    └── C1_Mask_Movie_Junctions.txt
```

---

## Parameters Reference

### Time Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `time_step` | float | 0.05 | Simulation timestep (seconds) |
| `total_time` | float | 10.0 | Total simulation duration (seconds) |
| `seed` | int | 42 | Random seed for reproducibility |

### Geometry Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `spindle_half_length` | float | 0.4 | Half-length of spindle (μm) |
| `spindle_width` | float | 0.32 | Width of spindle (μm) |
| `spindle_angle_init` | float | 0.0 | Initial spindle angle (radians) |
| `cell_radius_a` | float | varies | Cell semi-major axis (μm) |
| `cell_radius_b` | float | varies | Cell semi-minor axis (μm) |

### Microtubule Parameters

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `n_astral_mts` | int | 100 | 50-200 | Number of MTs per pole |
| `mt_mean_length` | float | 1.0 | 0.5-10 | Mean MT length (μm) |
| `mt_length_stdev` | float | 1.0 | 0.1-2 | Length variability |
| `length_distribution` | str | "gamma" | - | "gamma" or "constant" |
| `initial_state_distribution` | str | "random" | - | "random" or "growing" |

### Force Parameters

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `pull_force` | float | 5.0 | 0-20 | Motor pulling force (pN) |
| `push_force` | float | 0.0 | 0-20 | MT pushing force (pN) |
| `repulsion_force` | float | 0.0 | 0-10 | Cortex repulsion (pN) |

### Dynamics Parameters

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `growth_rate` | float | 0.013 | 0.005-0.03 | MT polymerization (μm/s) |
| `shrink_rate` | float | -0.027 | -0.05 to -0.01 | MT depolymerization (μm/s) |
| `catastrophe_rate` | float | 0.021 | 0.005-0.1 | Growth→shrink (1/s) |
| `rescue_rate` | float | 0.029 | 0.005-0.1 | Shrink→growth (1/s) |
| `bind_rate` | float | 0.03 | 0.01-0.1 | Motor binding (1/s) |
| `unbind_rate` | float | 0.02 | 0.01-0.1 | Motor unbinding (1/s) |

### Physical Constants

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `viscosity` | float | 100.0 | Cytoplasmic viscosity (Pa·s) |
| `bending_rigidity` | float | 2.0 | MT bending rigidity (pN·μm²) |
| `min_cortex_distance` | float | 0.01 | Minimum spindle-cortex gap (μm) |
| `max_interaction_distance` | float | 0.02 | Motor binding distance (μm) |
| `mt_min_length` | float | 0.01 | Minimum MT length before nucleation (μm) |

### Force Generator Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fg_density` | int | 100 | Number of cortical binding sites |
| `fg_density_basal` | int | 100 | Basal FGs (neuroblast only) |
| `fg_density_apical` | int | 100 | Apical FGs (neuroblast only) |
| `fg_density_anterior` | int | 40 | Anterior FGs (C. elegans) |
| `fg_density_posterior` | int | 60 | Posterior FGs (C. elegans) |

### Output Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_dir` | Path | Path("./output") | Output directory |
| `save_plots` | bool | True | Save visualization frames |
| `plot_interval` | int | 1 | Save every N steps |
| `show_force_vectors` | bool | False | Show force arrows in plots |

---

## Running Simulations

### Single Simulation

```python
from spindle import CElegansSpindleConfig, simulate

config = CElegansSpindleConfig(
    time_step=0.05,
    total_time=10.0,
    pull_force=5.0
)

results = simulate(config, run_id=1)
```

### Multiple Replicates

```python
# Run 10 replicates with different random seeds
for run_id in range(1, 11):
    config = CElegansSpindleConfig(
        time_step=0.05,
        total_time=10.0,
        seed=run_id,  # Different seed each time
        save_plots=(run_id == 1)  # Only save plots for first run
    )
    
    results = simulate(config, run_id=run_id)
    print(f"Run {run_id} complete: final angle = {results['angles'][-1]:.1f}°")
```

### Parameter Sweep

```python
# Test different pulling forces
pull_forces = [0.0, 2.5, 5.0, 7.5, 10.0]

for pull in pull_forces:
    config = CElegansSpindleConfig(
        time_step=0.05,
        total_time=10.0,
        pull_force=pull,
        save_plots=False
    )
    
    results = simulate(config, run_id=1)
    print(f"Pull={pull} pN: final angle={results['angles'][-1]:.1f}°")
```

### Using SLURM Array Jobs

For HPC clusters with SLURM:

**submit_job.sh:**
```bash
#!/bin/bash
#SBATCH --array=1-50
#SBATCH --time=2:00:00
#SBATCH --mem=4G

python run_simulation.py config.json
```

**run_simulation.py:**
```python
import os
import sys
from pathlib import Path
from spindle import load_config, simulate

# Get run ID from SLURM
run_id = int(os.getenv('SLURM_ARRAY_TASK_ID', 1))

# Load config
config_file = Path(sys.argv[1])
config = load_config(config_file)

# Update seed for this run
config.seed = run_id

# Run simulation
results = simulate(config, run_id=run_id)
```

Submit:
```bash
sbatch submit_job.sh
```

### Saving and Loading Configurations

**Save configuration:**
```python
from spindle import CElegansSpindleConfig, save_config

config = CElegansSpindleConfig(
    time_step=0.05,
    total_time=10.0,
    pull_force=5.0
)

save_config(config, "my_config.json")
```

**Load and run:**
```python
from spindle import load_config, simulate

config = load_config("my_config.json")
results = simulate(config, run_id=1)
```

**Command line:**
```bash
python spindle.py my_config.json
```

---

## Output and Analysis

### Output Directory Structure

```
output/
└── SM_{cell_type}_run_{run_id}/
    ├── config.json                      # Configuration used
    ├── {cell_type}_0.pdf               # Initial state
    ├── {cell_type}_1.pdf               # Frame 1
    ├── {cell_type}_2.pdf               # Frame 2
    ├── ...
    └── {cell_type}_data.xlsx           # Results data
```

### Results Dictionary

The `simulate()` function returns a dictionary:

```python
results = {
    'config': config,              # Configuration object
    'trajectory': [state1, ...],   # List of SimulationState objects
    'angles': np.array([...]),     # Spindle angles (degrees)
    'centers': np.array([[...]]),  # Center of mass positions (μm)
    'final_state': state_final,    # Final state object
    'push_forces': [forces1, ...], # Push forces each timestep
    'pull_forces': [forces1, ...]  # Pull forces each timestep
}
```

### SimulationState Object

Each state in the trajectory contains:

```python
state.time                  # Current time (s)
state.cell_boundary         # Cell outline (N, 2)
state.fgs                   # Force generator positions (M, 2)
state.spindle_poles         # Pole positions (2, 2)
state.spindle_angle         # Angle (radians)
state.astral_mts            # MT coordinates (2, n_mts, discr, 2)
state.astral_angles         # MT angles (2, n_mts)
state.mt_states             # Growth states: +1=growing, -1=shrinking (2, n_mts)
state.mt_lengths            # MT lengths (2, n_mts)
state.which_bind            # Binding status: 1=bound, 0=free (2, n_mts)
state.which_push            # Pushing status: 1=pushing, 0=not (2, n_mts)
state.free_fgs              # FG availability (M,)
state.astral_which_fg       # MT→FG assignments (M, 2)
state.velocity_com          # Velocity of center of mass (2,)
```

### Excel Output

The Excel file contains 8 sheets:

1. **Angle** - Spindle angle over time (degrees)
2. **Center** - Distance of spindle COM from origin (μm)
3. **Y-pos** - Y-coordinate of spindle COM (μm)
4. **N_pull** - Number of pulling (bound) MTs
5. **N_push** - Number of pushing MTs
6. **Ratio** - Pulling force ratio: 100 × F_pull / (F_pull + F_push)
7. **Pull_t** - Total pulling force magnitude (pN)
8. **Push_t** - Total pushing force magnitude (pN)

Each sheet has column `Run X` where X is the run_id.

### Reading Results

```python
import pandas as pd
from pathlib import Path

# Read Excel file
excel_file = Path("output/SM_celegans_spindle_run_1/celegans_spindle_data.xlsx")
data = pd.read_excel(excel_file, sheet_name=None, engine='openpyxl')

# Access individual sheets
angles = data['Angle']['Run 1'].values
n_pull = data['N_pull']['Run 1'].values
n_push = data['N_push']['Run 1'].values
ratio = data['Ratio']['Run 1'].values

print(f"Final angle: {angles[-1]:.1f}°")
print(f"Mean pulling MTs: {n_pull.mean():.1f}")
print(f"Mean force ratio: {ratio.mean():.1f}%")
```

### Analyzing Multiple Runs

```python
import pandas as pd
import numpy as np

# Combine data from multiple runs
all_angles = []

for run_id in range(1, 11):
    excel_file = Path(f"output/SM_celegans_spindle_run_{run_id}/celegans_spindle_data.xlsx")
    df = pd.read_excel(excel_file, sheet_name='Angle', engine='openpyxl')
    all_angles.append(df.iloc[:, 0].values)

all_angles = np.array(all_angles)

# Calculate statistics
mean_angle = np.mean(all_angles, axis=0)
std_angle = np.std(all_angles, axis=0)
final_angles = all_angles[:, -1]

print(f"Final angle: {final_angles.mean():.1f} ± {final_angles.std():.1f}°")
```

### Visualization Analysis

```python
import matplotlib.pyplot as plt

# Plot angle trajectory
times = np.arange(len(angles)) * 0.05  # Assuming 0.05s timestep

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# Angle over time
ax1.plot(times, angles)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Spindle Angle (°)')
ax1.set_title('Spindle Orientation')

# MT counts
ax2.plot(times, n_pull, label='Pulling')
ax2.plot(times, n_push, label='Pushing')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Number of MTs')
ax2.set_title('Microtubule Forces')
ax2.legend()

plt.tight_layout()
plt.savefig('analysis.pdf')
```

---

## Advanced Usage

### Custom Force Generator Distribution

For custom FG placement, modify the `make_fgs()` function or create your own:

```python
import numpy as np
from spindle import CElegansSpindleConfig, simulate

# Create custom FG distribution
def custom_fgs():
    # Example: FGs only on right half of cell
    angles = np.linspace(-np.pi/2, np.pi/2, 50)
    a, b = 2.5, 1.5
    fgs = np.column_stack([
        a * np.cos(angles),
        b * np.sin(angles)
    ])
    return fgs

# You would need to modify the simulation code to use this
# This is just to illustrate the concept
```

### Modifying MT Dynamics

To change how MTs behave, adjust the dynamics parameters:

```python
# Very stable MTs (low catastrophe, high rescue)
config = CElegansSpindleConfig(
    catastrophe_rate=0.005,  # Rare catastrophe
    rescue_rate=0.1,          # Frequent rescue
    growth_rate=0.02,         # Fast growth
    shrink_rate=-0.01         # Slow shrinking
)

# Highly dynamic MTs
config = CElegansSpindleConfig(
    catastrophe_rate=0.1,     # Frequent catastrophe
    rescue_rate=0.005,        # Rare rescue
    growth_rate=0.01,         # Slow growth
    shrink_rate=-0.05         # Fast shrinking
)
```

### Testing Force Balance

```python
# Pulling-dominated
config_pull = CElegansSpindleConfig(
    pull_force=10.0,
    push_force=2.0
)

# Pushing-dominated
config_push = CElegansSpindleConfig(
    pull_force=2.0,
    push_force=10.0
)

# Balanced
config_balanced = CElegansSpindleConfig(
    pull_force=5.0,
    push_force=5.0
)

# Compare results
results_pull = simulate(config_pull, run_id=1)
results_push = simulate(config_push, run_id=2)
results_balanced = simulate(config_balanced, run_id=3)
```

### Changing Cell Geometry

For elliptical cells, adjust aspect ratio:

```python
# Elongated cell
config = FollicleEpithelialConfig(
    cell_radius_a=1.0,   # Long axis
    cell_radius_b=0.5    # Short axis
)

# Circular cell
config = FollicleEpithelialConfig(
    cell_radius_a=0.75,
    cell_radius_b=0.75
)
```

For C. elegans, adjust superellipse parameter:

```python
# More rectangular
config = CElegansSpindleConfig(
    superellipse_n=1.5
)

# More diamond-shaped
config = CElegansSpindleConfig(
    superellipse_n=3.0
)
```

### Exporting for Further Analysis

```python
import pickle

# Save full results
with open('results.pkl', 'wb') as f:
    pickle.dump(results, f)

# Load later
with open('results.pkl', 'rb') as f:
    results = pickle.load(f)

# Extract specific data
trajectory = results['trajectory']
final_state = results['final_state']

# Save trajectory as NumPy arrays
np.save('spindle_poles.npy', 
        np.array([s.spindle_poles for s in trajectory]))
np.save('mt_lengths.npy',
        np.array([s.mt_lengths for s in trajectory]))
```

---

## Troubleshooting

### Common Issues

**1. Spindle moves outside cell**

**Symptoms:** Warning messages about boundary violations, spindle disappears

**Solutions:**
- Increase `min_cortex_distance` (default 0.01 μm)
- Reduce `max_step` (default 0.1)
- Lower force magnitudes
- Check initial spindle position with `spindle_angle_init`

```python
config = CElegansSpindleConfig(
    min_cortex_distance=0.02,  # Larger buffer
    max_step=0.05,              # Smaller max step
    pull_force=3.0              # Weaker forces
)
```

**2. Simulation very slow**

**Symptoms:** Takes many minutes for short simulations

**Solutions:**
- Reduce `n_astral_mts`
- Increase `time_step` (but check accuracy)
- Reduce `total_time`
- Set `save_plots=False` or increase `plot_interval`

```python
config = CElegansSpindleConfig(
    n_astral_mts=50,      # Fewer MTs
    time_step=0.1,         # Larger timestep
    plot_interval=10,      # Save fewer frames
    save_plots=False       # Skip plotting
)
```

**3. MTs not binding**

**Symptoms:** `which_bind` always 0, no pulling forces

**Solutions:**
- Check `pull_force > 0`
- Increase `fg_density`
- Increase `max_interaction_distance`
- Check `bind_rate > 0`

```python
config = CElegansSpindleConfig(
    pull_force=5.0,                  # Must be non-zero
    fg_density=150,                  # More FGs
    max_interaction_distance=0.03,   # Larger binding radius
    bind_rate=0.05                   # Faster binding
)
```

**4. Results don't match expected behavior**

**Solutions:**
- Check random seed for reproducibility
- Verify parameter units (μm, seconds, pN)
- Compare with example configurations
- Check force balance with `show_force_vectors=True`

```python
config = CElegansSpindleConfig(
    seed=42,                     # Fixed seed
    show_force_vectors=True      # Visualize forces
)
```

**5. Zebrafish cells: FileNotFoundError**

**Symptoms:** Cannot find cell images or Movie_info.xlsx

**Solutions:**
- Check data paths are correct
- Verify `endo_index` is valid (1-15)
- Ensure images exist for selected cell
- Check file naming convention matches

```python
from pathlib import Path

config = ZebrafishEndoConfig(
    endo_index=1,
    movie_info_file=Path("./data/Movie_info.xlsx"),
    cell_images_dir=Path("./data/cells"),
    spindle_images_dir=Path("./data/spindles")
)

# Verify files exist
print(config.movie_info_file.exists())
print(config.cell_images_dir.exists())
```

**6. MT lengths negative or unrealistic**

**Solutions:**
- Check `mt_mean_length` and `mt_length_stdev` are positive
- Verify `growth_rate` and `shrink_rate` have correct signs
- Check `mt_min_length` threshold

```python
config = CElegansSpindleConfig(
    mt_mean_length=1.0,      # Must be positive
    mt_length_stdev=0.3,     # Must be positive
    growth_rate=0.013,       # Must be positive
    shrink_rate=-0.027,      # Must be negative
    mt_min_length=0.01       # Minimum before renucleation
)
```

### Debugging Tips

**Enable verbose output:**
```python
# The simulate() function prints progress every 100 steps
# Watch for warnings about boundary conditions
results = simulate(config, run_id=1)
```

**Check intermediate states:**
```python
results = simulate(config, run_id=1)

# Examine state at different times
initial = results['trajectory'][0]
midpoint = results['trajectory'][len(results['trajectory'])//2]
final = results['final_state']

print(f"Initial angle: {np.rad2deg(initial.spindle_angle):.1f}°")
print(f"Midpoint angle: {np.rad2deg(midpoint.spindle_angle):.1f}°")
print(f"Final angle: {np.rad2deg(final.spindle_angle):.1f}°")

print(f"Initial pulling MTs: {np.sum(initial.which_bind)}")
print(f"Final pulling MTs: {np.sum(final.which_bind)}")
```

**Visualize forces:**
```python
config = CElegansSpindleConfig(
    show_force_vectors=True,  # Show all force arrows
    plot_interval=1            # Save every frame
)
results = simulate(config, run_id=1)
# Inspect PDF files to see force vectors
```

**Test with minimal configuration:**
```python
# Simplest possible simulation
config = CElegansSpindleConfig(
    time_step=0.1,         # Large timestep
    total_time=5.0,        # Short simulation
    n_astral_mts=20,       # Few MTs
    fg_density=30,         # Few FGs
    save_plots=True,
    plot_interval=1
)
results = simulate(config, run_id=1)
```

### Getting Help

If issues persist:

1. Check parameter ranges in [Parameters Reference](#parameters-reference)
2. Compare with working examples
3. Verify input data format (for zebrafish)
4. Check that all dependencies are up to date
5. Try with default parameters first

---

## Best Practices

### Reproducibility

Always set seeds for reproducible results:

```python
config = CElegansSpindleConfig(
    seed=42,
    # ... other parameters
)
```

### Performance

For parameter sweeps, disable plotting:

```python
config = CElegansSpindleConfig(
    save_plots=False,  # Much faster
    # ... other parameters
)
```

### Data Management

Organize outputs by experiment:

```python
from pathlib import Path

config = CElegansSpindleConfig(
    output_dir=Path("./results/experiment1/sweep1"),
    # ... other parameters
)
```

### Documentation

Save configuration with results:

```python
from spindle import save_config

save_config(config, results_dir / "config_used.json")
```



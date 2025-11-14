# Spindle Positioning Simulation - User Guide

A computational model for simulating mitotic spindle positioning across multiple cell types during cell division.

## Quick Start

### Single Simulation
```bash
# Basic follicle epithelial cell simulation
python spindle.py --cell_type follicle_epithelial

# With custom parameters
python spindle.py --cell_type follicle_epithelial --n_astral_mts 100 --pull_force 8.0
```

### In case you want run a bunch of simulations with SLURM
```bash
# Edit parameter arrays in multi_param_submit.sh, then run:
./multi_param_submit.sh
```
# Adjust SLURM partition, time etc in slurm_submit.sh

## Cell Types & Key Parameters

### Follicle Epithelial (`follicle_epithelial`)
*Drosophila* follicular epithelium with basal/lateral force generators
```bash
python spindle.py --cell_type follicle_epithelial \
  --n_astral_mts 100 \           # Microtubules per pole
  --fg_density 100 \             # Number of force generators
  --pull_force 5.0 \             # Pulling force (pN)
  --spindle_half_length 0.45      # Half-length of spindle
```

### Neuroblast (`neuroblast`)
*Drosophila* neuroblasts with apical force generator enrichment
```bash
python spindle.py --cell_type neuroblast \
  --n_astral_mts 100 \
  --fg_density_apical 100 \      # Apical force generators
  --fg_density_basal 100 \       # Basal force generators  
  --pull_force 5.0
```

### C. elegans PNC (`celegans_pnc`)
Posterior-enriched force generators, superellipse geometry
```bash
python spindle.py --cell_type celegans_pnc \
  --n_astral_mts 100 \
  --fg_density_anterior 40 \     # Anterior force generators
  --fg_density_posterior 60 \    # Posterior force generators
  --mt_mean_length 9.0           # Longer MTs than other systems
```

### C. elegans Spindle (`celegans_spindle`)
Uniform force distribution, superellipse geometry
```bash
python spindle.py --cell_type celegans_spindle \
  --n_astral_mts 100 \
  --fg_density_anterior 40 \
  --fg_density_posterior 60 \
  --push_force 5.0               # Include pushing forces
```

### Zebrafish Endothelial (`zebrafish_endo`)
Live imaging data with dynamic cell shapes
```bash
python spindle.py --cell_type zebrafish_endo \
  --endo_index 11 \              # Cell index from experimental data
  --fg_density 10 \              # Density per unit perimeter (not count!)
  --fg_distribution uniform      # or "junctions" (requires junction data)
```

## Important Parameter Notes

### Microtubules (`n_astral_mts`)
- **Per pole** - if you set 100, total MTs = 200 (100 per spindle pole)
- Typical ranges: 50-200 per pole
- More MTs = stronger forces but slower simulation

### Force Generators (`fg_density`)
- **Zebrafish**: Actual density (FGs per unit perimeter) - default 10
- **All others**: Actual count of force generators - default 100


### Forces
- `pull_force`: Magnitude when MT binds to cortical motor (pN)
- `push_force`: Magnitude when MT pushes against cortex (pN)  
- Typical: 0-10 pN range

### Timing
- `time_step`: Integration timestep (s) - default 0.05
- `total_time`: Simulation duration (s) 
- Smaller timesteps = more accurate but slower


## Testing a range of values for selected parameters

### Single Parameter
Edit `multi_param_submit.sh`:
```bash
# Change these arrays:
n_astral_mts=(50 100 200)        # Parameter to sweep
pull_force=(5)                   # Keep others fixed
push_force=(0.0)
spindle_half_length=(0.45)
```

### Multi-Parameter
```bash
# Test different combinations:
n_astral_mts=(100 200)
pull_force=(3.0 5.0 8.0)
push_force=(0.0 2.0)
spindle_half_length=(0.4 0.5)
# Results in 2×3×2×2 = 24 parameter combinations
```

### Custom Parameters
Add any parameter to the `args` line:
```bash
args="--cell_type $CELL_TYPE --n_astral_mts $n_mts --viscosity 150 --time_step 0.01"
```

## Output Structure

### Single Run
```
output/SM_follicle_epithelial_run_1/
├── config.json                    # Parameters used
├── follicle_epithelial_*.pdf      # Simulation snapshots
└── follicle_epithelial_data.xlsx  # Time series data
```

### Parameter Sweep
```
output/2024-11-12/                 # Date of run
└── follicle_epithelial_MT_100_pull_5_push_0.0_SL_0.45_test_1/
    ├── config.json
    ├── *_images/
    │   ├── *_task_1/              # Run 1 plots
    │   └── *_task_20/             # Run 20 plots
    └── *_stats/
        ├── task_1.xlsx            # Run 1 data
        └── task_20.xlsx           # Run 20 data
```

## Excel Data Format

Each Excel file contains sheets:
- **Angle**: Spindle angle over time (degrees)
- **N_pull**: Number of pulling microtubules
- **N_push**: Number of pushing microtubules  
- **Pull_t**: Total pulling force (pN)
- **Push_t**: Total pushing force (pN)
- **Ratio**: Pull/(Pull+Push) percentage

## SLURM Management

```bash
# Submit parameter sweep
./multi_param_submit.sh

# Check jobs
squeue -u $USER

# Cancel all jobs
scancel -u $USER

```

# Spindle Orientation and Positioning Model

A computational model for simulating mitotic spindle orientation and positioning during cell division. 

## Installation

```bash
pip install numpy pandas openpyxl matplotlib opencv-python shapely scipy
```

## Quick Start

```python
from spindle import CElegansSpindleConfig, simulate

# Create and run simulation
config = CElegansSpindleConfig(
    time_step=0.05,
    total_time=10.0,
    n_astral_mts=100,
    pull_force=5.0,
    push_force=5.0
)

results = simulate(config, run_id=1)
```

Or use built-in examples:

```python
from spindle import example_celegans_spindle
results = example_celegans_spindle()
```

## Supported Cell Types

- **Follicular Epithelial** (`FollicleEpithelialConfig`) - *Drosophila* egg chamber
- **Neuroblast** (`NeuroblastConfig`) - *Drosophila* with apical force generators
- **C. elegans PNC** (`CElegansPNCConfig`) - Posterior-enriched motors
- **C. elegans Spindle** (`CElegansSpindleConfig`) - Uniform distribution
- **Zebrafish Endothelial** (`ZebrafishEndoConfig`) - Live imaging data

## Output

Results are saved to `output/SM_{cell_type}_run_{id}/`:
- `config.json` - Simulation parameters
- `{cell_type}_{step}.pdf` - Visualization frames
- `{cell_type}_data.xlsx` - Time series data (angle, forces, MT counts)

# Spindle Orientation Model

A simulation of mitotic spindle orientation inside a cell cortex, driven by astral
microtubules (pushing) and cortical force-generator motors (pulling). Supports three
cell/organism systems: fly follicular epithelium/neuroblast, *C. elegans* (pronuclei
centering and metaphase/anaphase spindle positioning), and zebrafish embryonic cells
(using real tracked-movie cell shapes).

## Contents

- `spindle_model.py` — the model (single file, no build step).
- `Spindle_Simulator.ipynb` — widget-driven notebook for running a single simulation
  interactively (works locally or on Google Colab).
- `data/` — tracked-movie inputs for the zebrafish cell type (cell shapes, spindle
  masks, junction coordinates, movie metadata).

## Requirements

Python 3 with `numpy`, `matplotlib`, `opencv-python`, `pandas`, `openpyxl`, `shapely`,
`scipy`, and `ipywidgets` (only needed for the notebook UI).

## Quickstart (notebook)

Open `Spindle_Simulator.ipynb` and run the cells top to bottom — pick a cell type,
adjust the settings, click **Run simulation**.

## Quickstart (Python)

```python
from spindle_model import run_simulation

results = run_simulation(
    cell_type="FE",       # "FE" | "celegans" | "endo"
    motor_density=10,     # motor count (FE/celegans) or density (endo)
    n_astral_mt=50,
    push=1, pull=0,
    spindle_length=0.8,   # model units; FE and celegans "spindle" mode only
    total_time=600,       # seconds; None -> cell-type default
)
```

Key kwargs: `neuroblast=True` (FE apical-only placement), `celegans_mode="PNC"|"spindle"`,
`endo_cell_number=<n>` (required for `cell_type="endo"`; see `list_available_endo_cells()`),
`initial_angle`, `time_step`, `seed`, `advanced={...}` (override physical constants). Returns
a dict of result arrays (angle, position, push/pull counts vs. time) and writes an `.xlsx` to
`data_folder_path`.

## Quickstart (command line / Slurm)

```
python spindle_model.py <test_folder_path> <data_folder_path> <run_name>
```

`run_name` encodes the cell type (via a `cell_<Type>` token) and every parameter as regex
tokens, e.g.:

```
spindle_9214_cell_FE_MUD_10_MT_50_push_1_pull_0_ts_0.05_AL_1_state_default_angle_45_SL_1.8
```

`MUD`=motor density, `MT`=astral MT count, `push`/`pull`=on-off flags, `ts`=time step,
`AL`=astral MT length distribution, `state`=MT state-init mode, `angle`=initial spindle
angle, `SL`=spindle length. For `celegans` runs, the name must also contain `_PNC_` or
`_position_`; for zebrafish/`endo` runs, a `cell_<n>` token selects the tracked cell. Reads
`SLURM_ARRAY_TASK_ID`/`SLURM_JOB_ID` if set, otherwise runs as a plain local command.

## Cell types

- **FE** — fly follicular epithelium / neuroblast (idealized elliptical cell).
- **celegans** — *C. elegans* embryo; `PNC` (pronuclei centering) or `spindle`
  (metaphase/anaphase positioning) sub-modes.
- **endo** — zebrafish embryonic cell, using a real cell shape tracked from a movie
  (requires the matching entry under `data/`).

Units: lengths are in model units (1 model unit = 10 µm); time is in seconds.

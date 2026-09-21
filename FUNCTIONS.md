# spindle_model.py — Function Reference

This is a function-by-function reference for `spindle_model.py`, a single-file simulation of
mitotic spindle orientation inside a cell cortex.

## Global-variable coupling (read this first)

Many functions below do **not** take all their inputs as parameters — they read module-level
globals such as `cell_type`, `a`, `b`, `r`, `w`, `spread`, `config`, `time_step`, `total_time`,
`prob_dyn_bind`, `prob_dyn_unbind`, `prob_catastr`, `prob_rescue`, `AL`, `scale`, `NB`, `model`,
`c` (endo tracked-cell number), etc. These are set as a batch by `run_simulation(...)` (via
`globals()[...] = value`, since a function can't `global` a name that's also one of its own
parameters) before the relevant `simulate_*`/`plot_*` functions are called. Practically this
means:

- Most functions here are **not pure** — the same call can behave differently (or crash with a
  `NameError`) depending on what ran before it in the same process.
- `SimulationParameters()` must be instantiated only after `cell_type` (and, for `FE`/`endo`,
  `RC_rate`/`GS_rate`) are already set — several of its fields (astral MT length distribution,
  catastrophe/rescue/growth/shrink rates) branch on those globals.
- If you extract a function out of this file to reuse elsewhere, check whether it reads any of
  these globals before assuming it works standalone from its signature alone.

---

## 1. Shared / cell-type-agnostic helpers

Geometry, distance/perimeter math, and small stateless utilities used across all three cell
types.

### `angle_between_vectors(v1, v2)`
Angle in degrees between two 2D vectors, via `arccos` of the normalized dot product.

### `bounded_normal_random(mean, stdev)`
Draws from `random.gauss(mean, stdev)` repeatedly until the sample falls in `[0, 1]`. Used
historically for bounded random MT-length scaling (mostly superseded by the gamma-distributed
`orig_length` in `make_astral_MTs`).

### `calculate_partial_perimeter(coordinates, i, j)` / `calculate_perimeter(coordinates)`
Arc length of a closed polygon `coordinates` (N×2) between indices `i` and `j` inclusive (partial
perimeter), or around the whole loop (full perimeter). Used to convert a zebrafish FG "density"
parameter into an absolute point count along the real tracked-cell perimeter.

### `check_bind(i, j, astral, spindle_poles, spots, free_spots, astral_which_spot)`
If an astral MT tip (`astral`) is within `config.max_interact_dist` of the nearest free FG spot,
and a random draw succeeds against `prob_dyn_bind` (global), marks that spot occupied and returns
`bind=1`; otherwise `bind=0`. Mutates and returns `free_spots`/`astral_which_spot` bookkeeping
arrays. Superseded in the main loop by the combined `check_push_bind`/`check_push_bind_init`, but
still used in a couple of code paths.

### `check_push(a, b, astral, state, astral_angles, spindle_poles, cell)`
Returns `push=1` if an MT (assumed not already bound) has reached the cortex — i.e. its tip
distance from the pole matches the geometric intersection distance to `cell`, within floating
point tolerance — and is in the growing state (`state==1`).

### `check_push_bind_init(...)` / `check_push_bind(...)`
The combined bind-then-push check actually used during the simulation loop (`check_push_bind` is
the live one; `check_push_bind_init` is an earlier variant used at MT (re)birth). First checks
whether the MT tip is close enough to a free FG spot to bind (marking the spot taken and setting
`state=-1`); if not, checks whether it's close enough to the cortex to push (`state=1`). Reads
`config.max_interact_dist`/`config.push_dist` and (for `check_push_bind`) the global
`prob_dyn_bind`. `check_push_bind` additionally branches on `cell_type` for the "close enough to
cortex" test (`celegans` compares the vector difference; follicle/zebrafish compare the
difference of magnitudes).

### `check_push_junc(a, b, astral, bind, state, astral_angles, spindle_poles, cell, spots)`
Zebrafish-specific: suppresses pushing for MTs whose angle falls within the cell-cell junction
arc (`spots[0]`..`spots[-1]`), since junctions shouldn't generate cortical pushing.

### `check_spindle(spindle_poles, spindle_angle, cell, r, w)`
Validity check used throughout the force-balance step: builds the spindle body polygon
(`generate_spindle`) and returns `True` only if it stays at least `config.min_cortex_dist` inside
the cell polygon boundary (via Shapely `Polygon.buffer`/`.boundary.distance`).

### `close_contours(contours)`
Inserts evenly-spaced interpolated points between a contour's last and first point (spacing
≤ 0.005) so the loop closes without a long gap — used when preparing zebrafish tracked-cell
contours.

### `create_color_gradient(value)`
Maps `value` in `[0, 1]` to an RGBA color: green→white for `[0, 0.5]`, white→red for
`(0.5, 1]`. Used to color the pull/push-ratio indicator square in the per-frame plots.

### `create_simulation_directory(test_folder_path, name, task_id)`
Creates (if needed) and returns `(folder_name, new_dir_path)` for
`<test_folder_path>/<name>_<task_id>`, where per-frame PDFs get written.

### `delta_theta(a, b, theta)`
Re-parameterizes an ellipse's angular coordinate `theta` so that evenly-spaced *arc-length*
positions (rather than evenly-spaced angles) land correctly on an `a≠b` ellipse; used by
`make_fgs_follicle` when the cell is not circular.

### `distance_to_boundary(point, shape_coords)`
Minimum Euclidean distance from a point to any vertex of `shape_coords`. (Only used in a couple
of now-removed debug annotations; kept as a general utility.)

### `draw_vectors(ax, origin, forces, labels, colors, ...)`
Draws a set of labeled, arrow-normalized force vectors from a common `origin` on a Matplotlib
axis — used by the optional free-body-diagram overlay in the per-frame plots
(`config.show_vectors`).

### `enhance_cell(cell)`
Re-densifies a polygon's vertices so consecutive spacing never exceeds 0.005, by linearly
interpolating extra points into long edges. Used right after loading a raw tracked-cell contour.

### `find_number_after_cell(input_string)`
Regex-extracts the integer after `endo_` in a run name (the tracked zebrafish cell number) for
the CLI path.

### `find_tangential_and_normal_vectors(spindle_poles, spindle_angle, r, w, cell)`
Finds the closest point between the spindle body and the cell boundary, then returns the unit
tangent and (outward) normal vectors of the cell edge at that contact point. Used by
`slide_and_rotate_spindle`'s general approach (see below), though the live implementation now
searches along the raw velocity direction rather than strictly this tangent.

### `flip_angles(df_angle)` / `flip_angles_NB(df_angle)`
Post-processing on a full time series of spindle angles (in degrees): if the final (or, for the
NB variant, second) angle exceeds a threshold, reflects the whole series (`180 - angle` or
`180 + angle`) so all runs are reported on a consistent side, regardless of which way the spindle
happened to rotate.

### `gaussian_spaced_indices(x, y, n_samples, mean=None, std_dev=None)`
Generates `n_samples` indices between `x` and `y`, spaced according to a Gaussian quantile
function (denser near `mean`). Alternative to uniform spacing; not on the main FG-placement path
(that uses `gauss_points`/`get_sequence*` instead).

### `generate_spindle(spindle_poles, spindle_angle, r, w)`
Builds a 36-point 2D superellipse polygon representing the physical spindle body (semi-axes `r`,
`w`), centered at the midpoint of `spindle_poles` and rotated to `spindle_angle`. Central to all
cortex-containment checks (`check_spindle`, `project_spindle_to_feasible`).

### `get_astral_angle(astral)` / `get_spindle_angle(spindle)`
Angle (radians, normalized into `[0, 2π]`) of the vector from an MT's minus end to its plus end,
or between the two spindle poles, respectively.

### `get_contours(image_path)`
Loads a tracked-movie mask image (`cv2`), Gaussian-blurs and Canny-edge-detects it, extracts and
smooths the external contour, and returns it as an (N,2) array (y-flipped to match the model's
coordinate convention). Used for both cell and spindle masks in the zebrafish path.

### `get_sequence(n_points)` / `get_sequence_celegans(n_points)` / `get_sequence_top(n_points)` and `log_func*` helpers
A family of functions that generate non-uniformly-spaced angular sequences (via the paired
`log_func`/`log_func_celegans`/`log_func_top`/`log_func_top2` inverse-log mappings) for
distributing FGs non-linearly around a cortex. Largely superseded in the live FG-placement
functions by `gauss_points`, but still present/available.

### `grow_astralMT(a, b, angle, C, cell, orig_length)`
Advances one astral MT's plus end by one growth step (`config.growth_rate * time_step`) from pole
position `C` at `angle`, clamped to the cortex intersection if the naive grown length would
overshoot the cell boundary.

### `intersect_cell(a, b, angle, start_point, cell)`
Finds where a ray from `start_point` at `angle` first crosses the polygon `cell`, using Shapely
line/polygon intersection; picks the nearest intersection point if there are multiple, and falls
back to the geometric `intersect_cell_old` routine for degenerate intersection types. Returns
`(point, geom_type)`.

### `intersect_cell_old(a, b, angle, C, cell)`
A manual (non-Shapely) two-line-intersection fallback used when `intersect_cell`'s Shapely call
doesn't return a clean `Point`/`MultiPoint` result — walks `C` toward the cortex in `min(a,b)`
sized steps, then solves the two-nearest-vertices line intersection algebraically.

### `make_new_MT(i, j, a, b, MT, spindle_poles, astral_angles, cell, free_spots, astral_which_spot, orig_length)`
"Respawns" an astral MT that has shrunk below the minimum length: frees its old FG slot, resets
its length to 0, regrows it one step from the pole, and re-discretizes (`restructure`) the MT
point array.

### `normal_vector_to_ellipse(a, b, coordinate)`
Unit outward normal vector to the ellipse `x²/a² + y²/b² = 1` at a given `coordinate`, from the
gradient of the ellipse equation.

### `normalize_angles(angles)`
Wraps any scalar/array of angles (radians) into `[0, 2π)` via modulo.

### `point_in_polygon(point, polygon)`
Classic ray-casting point-in-polygon test; used to check whether a proposed spindle pole position
is still inside the cell.

### `restructure(beam)`
Re-discretizes a 2-point-or-more line segment (`beam`, first row = start, last row = end) into
`config.discr` evenly-spaced points via `linspace`. Used after any MT endpoint changes.

### `rotate_and_shift_parabola(x, y, theta, xi, yi)`
Rotates a set of `(x, y)` points by `theta` and translates by `(xi, yi)`. General-purpose
geometry helper.

### `sgn(x)`
Vectorized sign function (`np.where(x > 0, 1, -1)`) — note it returns `-1` for `x == 0` (no zero
case), used throughout the superellipse (celegans) shape math.

### `smoothen_cell(cell)`
Gaussian-smooths a contour's x/y coordinates independently (`scipy.ndimage.gaussian_filter1d`,
`sigma=10`) to remove tracked-movie pixel noise before resampling.

### `transfer_astro_which(new_length, b1)` / `transfer_free_spots(new_length, b1)`
Resize an `astral_which_spot`/`free_spots` bookkeeping array to `new_length`, copying over as much
of the old array `b1` as fits and zero-filling the rest. Used when the zebrafish FG spot count
changes as the tracked cell shape updates over time.

### `transform_junctions()`
Loads `data/junctions/C<c>_Mask_Movie_Junctions.txt` (global `c` = tracked cell number) and the
matching spindle mask image, rescales the junction coordinates into model units using the
spindle's own measured vs. tracked center offset, and returns the per-frame junction endpoint
pairs used by `make_fgs_zebrafish_junc`.

### `rotate_points(points, spindle_angle)`
Rotates an array of 2D points about their own centroid by `spindle_angle`.

### `slice_array(array, start, end, num_points)`
Returns `num_points` evenly-index-spaced rows from `array` between indices `start` and `end`
(inclusive) — used to pick a sub-arc of cell-boundary points for junction-restricted FG
placement.

### `distance_matrix(points1, points2)`
Pairwise Euclidean distance matrix between two point sets (N1×2, N2×2) → (N1×N2), via broadcasting
(no external dependency).

### `make_astral_MTs(params, cell, spindle_poles, spindle_angle, spots)`
Initializes the full astral-MT state for both spindle poles: angles (evenly spread across
`spread` around each pole's outward direction), initial lengths (gamma- or constant-distributed
depending on `config.length_MTs`), growth/shrink state (random or all-growing, depending on
`config.state_MTs`), then grows each MT one step and runs `check_push_bind` on it. Returns the
full `(astral_MTs, astral_angles, state, which_push, which_bind, free_spots, astral_which_spot,
orig_length, df_list2)` tuple threaded through the rest of the simulation. `params[0]` is
`[a, b, r, w]`, `params[1]` is the FG density/count, `params[2]` is the astral MT count per pole.

---

## 2. `SimulationParameters`

One config object holding every physical constant (rates, friction, viscosity, geometric
tolerances). **Must be constructed only after** the globals `cell_type` and (for `FE`/`endo`)
`RC_rate`/`GS_rate` are already set, since several fields branch on them.

### `SimulationParameters.__init__(self)`
Sets ~25 attributes: cortex vertex count (`number_of_sides`), astral MT length-distribution mode
(`length_MTs`, default `'gamma'` for all cell types), MT state-init mode (`state_MTs`), numeric
tolerances (`max_step`, `max_interact_dist`, `push_dist`, `min_cortex_dist`, `MT_min_length`,
`MT_max_length`, ...), then cell-type-branched microtubule dynamic-instability rates
(`catastr_rate`, `rescue_rate`, `shrink_rate`, `growth_rate` — `FE`/`endo` scale a fixed base rate
by the `RC_rate`/`GS_rate` multipliers; `celegans` hardcodes its own rates), motor
binding/unbinding rates, force magnitudes (`pull`, `push`, taken from the like-named globals),
rigidity/viscosity/friction (`EI`, `visc`, `mu_fric`), and a handful of feature toggles
(`pivoting`, `mobile_motors`, `force_velocity`, `show_vectors`, ...).

### `SimulationParameters.save_parameters_to_file(self, directory, filename, additional_params=None)`
Writes every attribute of `self` (plus an optional extra dict) as `key: value` lines to
`directory/filename`. No-ops (with a printed message) if the file already exists, so repeated
runs under the same name don't clobber the original parameter record.

---

## 3. Force/torque + spindle motion

The rigid-body force-balance step that moves and rotates the spindle each timestep, subject to
not penetrating the cortex.

### `gauss_points(N)`
Generates `N` values spaced via the inverse-CDF of a Gaussian (mean/sigma differ for `celegans`
vs. follicle, per a fit to Mud-protein distribution data), used as raw degree positions before
being filtered into an angular range by each `make_fgs_*` function.

### `chrom(x, y, angle)`
Returns two short line segments crossed like an "X" at `(x, y)`, oriented by `angle` — the
chromosome icon drawn in the per-frame plots. Segment half-length differs for `FE` vs. other cell
types.

### `new_spindle_poles(spindle_poles, spindle_angle, V, Omega)`
Advances the spindle poles/angle by one Euler step (`V * time_step`, `Omega * time_step`),
re-deriving the current half-length `r` from the pole separation for `celegans` (whose spindle can
elongate) or using the global `r` otherwise.

### `find_torque(spindle_poles, r, force1, force2)`
Computes net torque about the spindle's midpoint from the two poles' net forces, using the
transverse (perpendicular-to-spindle-axis) force components and their sign via a 2D cross
product. Re-derives `r` for `celegans` as in `new_spindle_poles`.

### `find_force(astral_MTs, spindle_poles, which_bind, which_push, v_c)`
For every astral MT, computes a pull vector (if bound to an FG — optionally velocity-dependent via
`config.force_velocity`/`config.v_0`) or a push vector (if pushing — zebrafish uses a
length-proportional force, follicle/celegans use an Euler-buckling-limited force
`EI·π²/length²`). Returns `(push_vecs, pull_vecs)`, each shaped `(2 poles, n_astro, 2)`.

### `check_push_bind(...)`
See section 1 (listed there since it's also used during MT initialization).

### `slide_and_rotate_spindle(spindle_poles, spindle_angle, V, Omega, cell, r, w, config, dt)`
Fallback used when the naive `new_spindle_poles` step would leave the spindle invalid
(`check_spindle` fails): binary-searches (10 iterations for zebrafish, 20 otherwise) for the
largest safe fraction of the intended translation `V*dt`, then — independently — the largest safe
fraction of the intended rotation `Omega*dt` at the translated position. Returns
`(new_poles, new_angle, moved: bool)`; if even a zero-fraction result is still invalid, returns the
original configuration unchanged. Prints verbose `[DEBUG]` traces of every bisection step.

---

## 4. Per-timestep update

### `update_astral_MTs(params, cell, spindle_poles, spindle_angle, delta_spindle_angle, astral_MTs, astral_angles, state, which_push, which_bind, spots, free_spots, astral_which_spot, orig_length, force_vector_1, force_vector_2, run)`
The per-MT state machine run every timestep for every astral MT on both poles. Refreshes each
MT's minus-end position to the (possibly moved) pole, then branches:
- **Bound MTs**: re-checks the cortex intersection and randomly unbinds (`prob_dyn_unbind`) or
  becomes geometrically infeasible (`MultiPoint` intersection); zebrafish additionally re-glues or
  re-finds the FG spot each refresh cycle depending on `config.mobile_motors`.
- **Unbound MTs**: rolls dynamic instability (`prob_catastr`/`prob_rescue` against the current
  `state`), then either shrinks (clamped to the cortex intersection, respawning via `make_new_MT`
  if it falls below `config.MT_min_length`) or grows (`grow_astralMT`), then re-runs
  `check_push_bind`.

Also builds a big per-MT `df_list2` log row (bind/push/state/force/length/angle/etc.) consumed by
`plot_simulate_*` for the optional raw-MT Excel export. Reads globals `spread`, `time_step`,
`prob_catastr`, `prob_rescue`, `prob_dyn_unbind`, `cell_type`, `config`.

### `move_spindle(params, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, cell, spots, which_push, which_bind, free_spots, astral_which_spot, orig_length, v_c, i)`
The main per-timestep driver: for zebrafish, periodically refreshes the tracked cell shape/FG
layout (`update_cell_zebrafish`, `make_fgs_zebrafish_*`) and resizes the bookkeeping arrays; computes
net push/pull force and torque (`find_force`, `find_torque`); derives translational/rotational
velocity via Stokes-drag-like formulas (`6πμr`, `8πμr³`); proposes new poles/angle
(`new_spindle_poles`), falling back to `slide_and_rotate_spindle` if invalid; then calls
`update_astral_MTs` for the new configuration. Returns the full updated state tuple plus the
pull/push force ratio. Prints verbose `[DEBUG]` traces.

### `move_severed_spindle(params, astral_MTs, astral_angles, state, spindle_poles, spindle_angle, cell, spots, which_push, which_bind, free_spots, astral_which_spot, orig_length, v_c, i)`
Used after a `celegans` anaphase "severance" event (spindle poles decoupled): moves each pole
independently as its own Stokes sphere (`move_sphere`, a nested helper implementing simple
drag-limited motion with a hard cortex-radius clamp), then still runs `update_astral_MTs` for both
poles.

---

## 5. Cell geometry over time (zebrafish only)

Zebrafish cortex shape is *read from tracked movies*, not simulated — these functions manage that.

### `resample_cell(traj, n_points=2400)`
Re-parameterizes a raw contour by arc length into exactly `n_points` evenly-spaced points
(`scipy.interpolate.interp1d`), then rotates the point order so index 0 is nearest the ray from
the origin at angle 0 (via `intersect_cell`) — giving a consistent starting vertex across frames.

### `get_real_cell(t_time=0)`
Loads the tracked cell mask image for the current movie frame (`data/cells/cell_<c>/Mask_<frame>.jpg`,
global `c`), extracts and rescales its contour (`get_contours`, `rescale[c-1]`), similarly loads
the matching spindle mask to compute the spindle half-length `r` and center offset, recenters the
cell contour on the spindle center, smooths (`smoothen_cell`), closes the loop
(`close_contours`), and resamples to 2400 points (`resample_cell`).

### `interpolate_contours(contour1, contour2, steps=5)`
Linearly interpolates `steps` intermediate contours between two equal-length point arrays — used
to produce sub-frame cell shapes between two tracked-movie frames.

### `update_cell_zebrafish(t_time)`
Gets the tracked cell shape at the current and next movie frame (`get_real_cell`), interpolates
`frame_rates[c-1]` steps between them, and returns the one matching the current sub-frame
position within that interval. Called periodically (once per simulated second) from
`move_spindle`.

### `project_spindle_to_feasible(spindle_poles, spindle_angle, cell, r, w, step=0.02, max_trans_only=0.1, rings=20, directions=16, rot_steps=10, rot_span_deg=20)`
Recovery routine for when the spindle has ended up outside the feasible region (can happen after
a zebrafish cell-shape refresh shrinks the cortex around it): Stage 1 searches nearby
translation-only offsets (rings of `step`-sized radii, capped at `max_trans_only`) for the nearest
feasible position at the same angle; Stage 2, if that fails, widens the search to combined
translation + rotation. Returns the original pose unchanged if nothing feasible is found within
the search budget. Verbose `[DEBUG]` prints throughout.

---

## 6. Cell-type-specific setup: FG placement and initial cell/spindle geometry

### FG (force-generator "motor") placement — `make_fgs_*(motors, a, b, cell, spindle_poles, frame)`
All variants take the same signature and return an (N,2) array of FG spot coordinates; they
differ only in where those spots go:
- **`make_fgs_follicle`** — a Gaussian-clustered band on the basolateral cortex (±60° from the
  poles), mirrored left/right; `motors` is a plain count.
- **`make_fgs_neuroblast`** — a single apical/top-only Gaussian-clustered band (30°–150°),
  replacing the basolateral placement (selected when the `NB` flag is set).
- **`make_fgs_celegans_pnc`** — FGs spread broadly (outside ±85°) around both anterior and
  posterior poles of the superellipse cortex, split 60/40 anterior/posterior.
- **`make_fgs_celegans_spindle`** — same idea but tighter cutoffs (outside ±90°) and a 40/60
  anterior/posterior split, matching the metaphase/anaphase posterior-shift geometry.
- **`make_fgs_zebrafish_junc`** — FGs restricted to the tracked cell-cell junction arc
  (`transform_junctions` + `calculate_partial_perimeter`); `motors` here is a **density**
  (multiplied by the junction arc length to get a point count).
- **`make_fgs_zebrafish_uf`** — FGs spread uniformly over the *entire* tracked perimeter;
  `motors` is again a density (multiplied by `calculate_perimeter(cell)`).

### Initial cell + spindle placement — `make_cell_follicle(params)` / `make_cell_celegans(params)` / `make_cell_zebrafish(params)`
All three build the initial cortex polygon, place the two spindle poles at the requested angle
(nudging the angle in `π/12` steps, up to 12 times, if the initial placement falls outside the
cell), call the matching `make_fgs_*` function(s) to place FGs, then call `make_astral_MTs` to
populate the astral MT arrays. They return the same 13-element tuple: `(cell, astral_MTs,
astral_angles, state, spindle_poles, spindle_angle, spots, which_push, which_bind, free_spots,
astral_which_spot, orig_length, df_list2)`.
- **`make_cell_follicle`** — cortex is a plain ellipse (`a`, `b`); dispatches to
  `make_fgs_neuroblast` or `make_fgs_follicle` depending on the `NB` flag.
- **`make_cell_celegans`** — cortex is a superellipse (fixed `a=2.5, b=1.5, n=2.2`); initial pole
  displacement differs slightly between `PNC` and `spindle` sub-modes; dispatches to
  `make_fgs_celegans_pnc`/`_spindle` based on the `model` global.
- **`make_cell_zebrafish`** — cortex comes from `get_real_cell()` (tracked movie); dispatches to
  `make_fgs_zebrafish_uf`/`_junc` based on the `model` global (`'ever'`/`'junc'`).

`params` throughout this file is the positional list `[[a, b, r, w], motor_density, n_astral_mt,
initial_spindle_angle, folder_name, output_dir_path]`.

---

## 7. Run-name parsing

Each function below regex-extracts one cell-type family's Slurm-style run-name token convention
into raw parameter values (mirrors what would have been 3 separate CLI scripts before the merge).
None of them read or set globals — they're pure string→tuple parsers.

### `extract_parameters_follicle(input_string)`
Extracts `MUD_<int>` (motor density), `MT_<int>` (astral MT count), `push_<float>`, `pull_<float>`,
`ts_<float>` (time step), `cell_<name>`, `AL_<letters>`/`AL_<float>` (length-distribution mode /
astral length), `state_<letters>`, `angle_<float>`, `SL_<float>` (spindle length), `mcd_<float>`,
and an `FE_NB` marker (neuroblast flag), plus `RC_rate_<float>`. Returns `(motor_density,
astral_MTs, push, pull, SL, time_step, cell_type, angle, NB, RC_rate)`. Any token not found in
the string yields `None` for that value (except `RC_rate`, which defaults to `1.0`).

### `extract_parameters_zebrafish(input_string)`
Same core tokens as above, plus a `junc`/`ever` mode marker and `GS_rate_<float>`. Returns
`(motor_density, astral_MTs, push, pull, SL, time_step, cell_type, angle, model, GS_rate, AL)`.

### `extract_parameters_celegans(input_string)`
Same core tokens; additionally requires the run name to contain `_PNC_` (→ `model="PNC"`,
`angle=90`) or `_position_` (→ `model="spindle"`, `angle=0`) — this is the actual mechanism the
`__main__` block relies on to pick the sub-model. Returns `(motor_density, astral_MTs, push,
pull, SL, time_step, cell_type, angle, model)`.

---

## 8. Per-frame plotting

Purely cosmetic — renders one Matplotlib PDF frame of the cell, spindle, color-coded astral MTs,
FGs, and (optionally) force vectors. Does not affect simulation results. All three take the same
long positional signature: `(cell, astral_MTs, astral_angles, state, spindle_poles, spindle_angle,
push_force, pull_force, spots, i, which_bind, which_push, v_c, n, ratio, params, borders)`.

### `plot_cell_follicle(...)` / `plot_cell_celegans(...)` / `plot_cell_zebrafish(...)`
Shared structure across all three:
- Draws the cortex outline, spindle poles (yellow dots), FG spots (salmon markers, labeled with
  astral MT count and either "Number of FGs" or — for zebrafish — "FG density"), and each astral
  MT color-coded by state: red = bound/pulling, dark green = pushing, cyan = shrinking, slate blue
  = free/growing.
- Legend entries: `Time = <n * time_step> s`, cortex angle, FG/force-vector labels, and net force
  per pole/total (`$F_{pole1}$`, `$F_{pole2}$`, `Total force`). Legend is forced fully opaque white
  (`facecolor='white', framealpha=1`) so it doesn't pick up a grey tint from the cortex/MT lines
  drawn underneath it.
- Optional free-body-diagram overlay of push/pull/repel/net force vectors when
  `config.show_vectors` matches a function-specific sentinel value.
- On-plot text: pull/push force ratio and pull/push MT counts, plus a color-graded square
  (`create_color_gradient`) indicating the pull/push balance.
- Draws chromosome icons (`chrom`) at the metaphase plate; `plot_cell_celegans` also draws a
  severance-event marker line once past `severance_time`.
- Saves to `<params[5]>/<params[4]>_<n+1>.pdf`.

Differences: `plot_cell_follicle` skips FG-count labeling when `pull == 0`; `plot_cell_celegans`
additionally draws the spindle body envelope for `model == 'PNC'` and shows spindle length in µm
(`20 * r`, since 1 model unit = 10 µm); `plot_cell_zebrafish` uses tracked-movie plot borders
(passed in via `borders`) instead of deriving them from the idealized cell, and labels the cortex
angle as "Spindle angle" relative to the tracked cell's long axis.

---

## 9. Simulation drivers

### `simulate_follicle(params)` / `simulate_celegans(params)` / `simulate_zebrafish(params)`
Each: builds the initial cell/spindle/MTs (`make_cell_*`), plots the initial frame if
`SAVE_FRAME_PLOTS`, then loops calling `move_spindle` (or, for `celegans` past
`severance_time` with `severance == 1`, `move_severed_spindle`) once per timestep until
`i > total_time / time_step`, accumulating a `plot_list` of per-step state snapshots and a
per-MT log `DataFrame`. Saves a frame plot every `save_interval_time` seconds
(2s for follicle/zebrafish, 1s for celegans), gated additionally by `task_id` (only the first
few Slurm array tasks save PDFs, to avoid flooding disk on large array jobs). Returns
`(plot_list, df)`.

### `plot_simulate_follicle(params)` / `plot_simulate_celegans(params)` / `plot_simulate_zebrafish(params)`
Thin wrappers: call the matching `simulate_*`, optionally write the raw per-MT log to
`<folder>.xlsx` (`task_id == 1` and `SAVE_RAW_MT_LOG`; zebrafish's variant doesn't currently write
this file), then extract `plot_list` into flat per-timestep series — `angle` (degrees,
`flip_angles`-normalized for follicle/celegans), `center`, `y-position`, push/pull MT counts,
pull/push ratio, and total pull/push force vectors. Returns an 8-tuple of these series, which
`run_simulation()` bundles into its result dict.

---

## 10. `run_simulation(...)`

The explicit-keyword-argument entry point (used by both the notebook and, indirectly, the CLI).
Validates `cell_type` (`"FE"`/`"celegans"`/`"endo"`) and that `endo_cell_number` is given for
`"endo"`; sets the shared globals (`SAVE_FRAME_PLOTS`, `SAVE_RAW_MT_LOG`, `task_id`, `cell_type`,
`push`, `pull`, `time_step`, `RC_rate`, `GS_rate`); seeds RNGs if `seed` is given; builds
`SimulationParameters()` and applies any `advanced={...}` attribute overrides; then branches by
`cell_type` to set every remaining cell-type-specific global (`a`, `b`, `r`, `w`, `spread`,
`total_time`, `AL`, `scale`, initial `sp_angle`, etc.) — including, for `celegans`/`FE`, resolving
`spindle_length` (`None` → 1.8 model units for `celegans` "spindle" mode, or 0.8 for `FE`) and, for
`FE`, raising `ValueError` if the resolved spindle length exceeds 0.95 model units (too long to
fit the `a=b=0.5` follicle cell). For `endo`, reads `cells_data.xlsx`/`Movie_info.xlsx` from
`DATA_DIR` to pull the tracked cell's initial angle, frame rate, and metaphase→anaphase duration.
Finally calls `create_simulation_directory`, saves the parameter file
(`SimulationParameters.save_parameters_to_file`), dispatches to the matching `plot_simulate_*`,
writes results to an `.xlsx` via `pd.ExcelWriter`, and returns a dict with `angle`, `count_push`,
`count_pull`, `center`, `ypos`, `ratio`, `pull_t`, `push_t`, `time_step`, `folder_name`, `run_dir`,
and `data_xlsx`.

---

## 11. CLI entry point (`if __name__ == "__main__":`)

Reads `test_folder_path`, `data_folder_path`, `run_name` from `sys.argv`; reads
`SLURM_ARRAY_TASK_ID`/`SLURM_JOB_ID` (defaulting to `1`/`None` outside an actual Slurm array job);
probes the cell type from a `cell_<Type>` token in `run_name`; dispatches to
`extract_parameters_follicle`/`_zebrafish`/`_celegans` to decode the rest of the tokens into
`run_kwargs`; and calls `run_simulation(...)` with those kwargs plus
`save_frame_plots=save_raw_mt_log=verbose=True` (batch runs always produce every artifact).

---

## 12. Notebook / UI helpers

Additive-only; not used by the CLI path. `build_ui()` imports `ipywidgets` lazily so importing
this module for batch use never requires it.

### `list_available_endo_cells()`
Returns the sorted list of tracked-cell numbers that can actually be run as `cell_type="endo"`:
present as `data/junctions/C<n>_Mask_Movie_Junctions.txt` **and** having a matching `cell_<n>`
sheet in `data/cells_data.xlsx`. Used to populate the notebook's tracked-cell dropdown with only
choices that will actually run.

### `plot_summary(results, title=None)`
Given a `run_simulation(...)` result dict, draws a 2×2 Matplotlib figure: spindle angle vs. time,
push/pull MT counts vs. time, spindle center x/y position vs. time, and the 2D center trajectory
(start/end markers). Returns the `Figure` without calling `plt.show()` (left to the caller/Jupyter
to display, to avoid blocking under non-inline backends).

### `build_ui()`
Builds and returns the full `ipywidgets` VBox used by `Spindle_Simulator.ipynb`: a cell-type
dropdown, initial-angle/spindle-length controls, fill-in number boxes for motor density and
astral MT count, a duration box, push/pull/save-output checkboxes, a seed control, run-label/
output-folder text boxes, an "Advanced" accordion for physical-constant overrides, and a Run
button. Internally defines:
- `_refresh_runtime_estimate()` — updates a rough wall-clock estimate from a per-cell-type
  `TIMING_FIT` (ms/step ≈ `base + slope * n_astral_mt`) fit to measured timings.
- `_refresh_endo_info()` — for `cell_type="endo"`, reads the tracked cell's initial angle and
  metaphase→anaphase duration from `Movie_info.xlsx` and pre-fills the duration box.
- `_refresh_fields(*_)` — swaps which controls are shown/labeled/ranged based on the selected
  cell type and sub-mode (e.g. spindle length in µm with an 8 µm default for `FE`, vs. model units
  with a 1.8 default for `celegans` "spindle" mode; motor field labeled "Motor density" only for
  `endo`).
- `_on_run_clicked(_)` — reads every widget's `.value`, converts FE's spindle length from µm back
  to model units (÷10), builds the `run_simulation(...)` kwargs for the selected cell type, runs
  it, and displays `plot_summary(...)` in an output area, or shows the exception message (e.g. the
  "too long" `ValueError` from `run_simulation`) on failure.

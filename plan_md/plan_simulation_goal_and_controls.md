# plan_simulation_goal_and_controls

## Goal of `2.18_simulation_fixed.R`
1. Run a simulation benchmark for sparse multivariate regression with sparse precision estimation.
2. Compare Prox-Newton (PN) and Prox-Gradient (PG) along a 60-point lambda path.
3. Evaluate convergence behavior under three stopping cases and export diagnostics, timings, and plots.

## Locked Protocol (Revised)
1. Outer tolerances:
   - `tol_PN_sub = 1e-7`
   - `tol_PG_sub = 1e-7`
2. PG minimum outer iterations:
   - `outer_pg$min_iter = 100`
   - PG stopping checks are enabled only when iteration index `m >= 100` (unless clipped by `max_iter`)
3. Inner tolerance ladder:
   - base `1e-8`
   - adaptive tightening by `/10`
   - floor `1e-10`
4. Three stop-mode cases:
   - `case_relative_change`: PN=`relative_change`, PG=`relative_change`
   - `case_oracle_change`: PN=`oracle_change`, PG=`oracle_change`
   - `case_norm_mode`: PN=`local_norm`, PG=`l2_step_norm`
5. Oracle source:
   - PN fixed-100 oracle path by lambda
   - oracle value per lambda is minimum loss across iterations `0..100`
   - same oracle vector is used by both PN and PG in `case_oracle_change`
6. Plot requirements:
   - 6 plots total (`3 cases x 2 metrics`)
   - metrics: training loss vs `log(lambda)`, iterations vs `log(lambda)`
   - PN/PG overlaid in each plot
   - each plot subtitle includes total path time for PN and PG
7. Timing export:
   - `solution_path_time_by_stop_mode_tol_1e-07.csv`
   - required columns: `case`, `method`, `total_elapsed_sec`, `tol_outer`, `tol_inner_base`, `tol_inner_floor`

## Task Mapping (Tasks 1–8)

### Task 1 + Task 8: Centralized controls and full visibility
1. Keep one config object: `SIM_CONFIG <- get_default_sim_config()`.
2. Keep all run-critical parameters in config blocks:
   - `simulation`, `lambda_path`, `outer_pg`, `outer_pn`
   - `pn_subproblem_admm`, `omega_prox_admm_pg`, `omega_prox_admm_pn`
   - `oracle_eval`, `numerics`, `benchmark`, `reporting`
3. Print active config each run with `print_sim_config(SIM_CONFIG)`.
4. Export parameter catalog each run:
   - `simulation_parameter_catalog.csv`
   - `simulation_parameter_catalog.rds`
5. Include numerical stability parameters explicitly (eigen floors, jitter, inverse/logdet eps, conditioning thresholds).

### Task 2: Oracle PN fixed-100 and non-decrease counting
1. Use `run_pn_oracle_fixed_iters(...)` over the full lambda path.
2. Run exactly 100 PN updates per lambda (101 recorded losses including iteration 0).
3. Oracle loss per lambda: `min(loss_hist[0:100])`.
4. Non-decrease criterion: `loss_new >= loss_prev - 1e-12`.
5. Export:
   - `oracle_pn_100iter_loss_summary.csv`
   - `oracle_pn_100iter_nondecrease_counts.csv`

### Task 3: PN and PG stop-mode controls
1. PN stop modes:
   - `local_norm`
   - `relative_change`
   - `oracle_change`
2. PG stop modes:
   - `relative_change`
   - `oracle_change`
   - `l2_step_norm` (`||d||_F`)
3. Stop logic is configurable (`any` / `all`), default `any`.

### Task 4: PN subproblem ADMM tolerance hierarchy and adaptive tightening
1. Use stricter inner tolerance than outer:
   - outer `1e-7`
   - PN subproblem ADMM base `1e-8`
2. If PN line-search fails to obtain descent after 3 halvings:
   - tighten active `tol_PN_ADMM <- max(tol_PN_ADMM/10, 1e-10)`
   - retry the same PN outer iteration immediately
3. Allow up to 3 refinement rounds per outer iteration.
4. If still no descent, stop without forced uphill acceptance.
5. Export refinement trace:
   - `pn_adaptive_tol_refinement_trace.csv`

### Task 5 and Task 7: Omega ADMM termination modes for PN/PG
1. Support two termination modes in `prox_psd_offdiag_l1(...)`:
   - `current`
   - `duality_gap`
2. Keep dedicated tolerances and max iterations for PN and PG Omega ADMM.
3. Export mode summaries:
   - `pn_omega_admm_mode_summary.csv`
   - `pg_omega_admm_mode_summary.csv`

### Task 6: PG outer update controls
1. PG outer tolerance and max iterations are config-driven.
2. `tol_PG_sub = 1e-7`.
3. `outer_pg$min_iter = 100` is enforced globally across all stop-mode cases.
4. PG stop checks are gated: for `m < min_iter_PG_outer`, no PG stop condition can terminate training.
5. If `min_iter_PG_outer > max_iter_PG_outer`, it is clipped to `max_iter_PG_outer`.
6. Case runner switches PG stop mode among:
   - `relative_change`
   - `oracle_change`
   - `l2_step_norm`

## Required Output Artifacts
1. Case plots (6 total):
   - `pn_pg_final_loss_vs_loglambda_case_relative_change_tol_1e-07.png`
   - `pn_pg_iters_vs_loglambda_case_relative_change_tol_1e-07.png`
   - `pn_pg_final_loss_vs_loglambda_case_oracle_change_tol_1e-07.png`
   - `pn_pg_iters_vs_loglambda_case_oracle_change_tol_1e-07.png`
   - `pn_pg_final_loss_vs_loglambda_case_norm_mode_tol_1e-07.png`
   - `pn_pg_iters_vs_loglambda_case_norm_mode_tol_1e-07.png`
2. Path-time export:
   - `solution_path_time_by_stop_mode_tol_1e-07.csv`
3. Keep existing diagnostics and run snapshots.

## Acceptance Criteria
1. Config snapshot shows `tol_PN_sub=1e-7` and `tol_PG_sub=1e-7`.
2. Inner tolerances initialize at `1e-8`, tighten by factor `10`, and never go below `1e-10`.
3. `run_config_snapshot_*` includes `min_iter_PG_outer` and shows value `100` (or clipped `max_iter` if `max_iter < 100`).
4. For each case, PG path iterations satisfy `outer_iters >= 100` for every lambda point when `max_iter >= 100`.
5. If `max_iter < 100`, PG uses `min_iter_PG_outer = max_iter` and stops only at cap or later checks.
6. PN and PG both run under all three case stop modes (method-appropriate norm modes).
7. `case_oracle_change` uses shared PN fixed-100 oracle vector for both methods.
8. Exactly 6 required plot files are produced with PN/PG overlays and timing subtitles.
9. `solution_path_time_by_stop_mode_tol_1e-07.csv` contains 6 rows (`3 cases x 2 methods`) with required columns.

## PN-first PG-target solver protocol (paper workflow)
1. Dedicated run profile:
   - `RUN_PROFILE=pn_first_pg_target`
2. Fixed setting for this workflow:
   - `case_norm_mode`
   - `tol_outer=1e-4`
   - `tol_inner_base=1e-5`
   - `tol_inner_floor=1e-8`
3. PN-first stage:
   - run PN across all 60 lambda paths
   - record PN final objective loss per path (`pn_final_loss`)
4. PG target-only stage:
   - run PG with strict target condition:
   - stop only when `PG loss < PN final loss - 1e-12`, or at hard cap
   - PG hard cap: `max_iter=10000`
   - PG stop modes are disabled in this profile
5. Time reporting for this profile:
   - per-path table with PN/PG elapsed time and tolerance metadata
   - total table with one row per method (`PN`, `PG`)
6. PN iteration distribution requirement:
   - frequency table by PN outer iterations (`pn_outer_iters`, `path_count`, `path_fraction`)
7. Required artifacts for this profile:
   - `pn_reference_final_loss_case_norm_mode_tol_1e-04.csv`
   - `pg_to_pn_target_path_case_norm_mode_tol_1e-04_max_iter_10000.csv`
   - `pn_pg_target_comparison_case_norm_mode_tol_1e-04_max_iter_10000.csv`
   - `pn_pg_time_per_path_case_norm_mode_tol_1e-04_max_iter_10000.csv`
   - `pn_pg_time_total_case_norm_mode_tol_1e-04_max_iter_10000.csv`
   - `pn_outer_iter_frequency_case_norm_mode_tol_1e-04.csv`
   - `run_config_snapshot_pn_first_pg_target_tol_1e-04.csv`

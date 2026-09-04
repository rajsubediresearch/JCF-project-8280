Run everything from the `forecasting_odemodels code/` folder, with that
folder on the MATLAB path.

Run the steps **in order** - each one reads files written by the
previous one.

# Part A - Practice run (Ebola scenario 1)

Nothing needs editing. All configuration is already set to $C = 20$
weeks, $H = 4$ weeks, $s = 4$ weeks, 5 windows.

### 1. Fit each model across all rolling windows

    batch_GLM_ebola_scenario1
    batch_GOM_ebola_scenario1
    batch_RICH_ebola_scenario1

Slowest step (bootstrap refits, `B = 300` per window per model). Writes
parameter estimates, forecast curves, bootstrap trajectories and
per-window performance files to `output/`.

### 2. Collect metrics

    collect_metrics_simple()

Writes `forecast_metrics_summary.mat` / `.csv` --- MAE, MSE, coverage,
WIS and AICc for every model-window pair. **This file supplies both JCF
inputs** (`WIS_calib` and `WIS_forecast`), so steps 1--2 never need
repeating for JCF.

### 3. Baseline ensemble (calibration-only weights)

    create_ensemble_simple()

Weights each model by $1/{WIS}^{cal}$, pools with the linear pool.
Writes `ensemble_results.mat` / `.csv` and
`output/Forecast-Ensemble-LinearPool-*.csv`.

### 4. Baseline plots

    plot_model_fit_simple()
    plot_forecasts_simple()

Figures go to `plots/`.

### 5. JCF ensembles

    run_jcf_sweep()

Runs $\lambda \in \{ 0,\ 0.25,\ 0.5,\ 0.75,\ \text{length},\ 1\}$, where
`'length'` means $\lambda = C/(C + H)$. Each run writes its own tagged
files (`ensemble_results_JCF-lam050.*`, `JCF-lamlen`, ...), so nothing
overwrites the baseline from step 3.

For a single $\lambda$ instead: `create_ensemble_jcf()` (edit
`config.lambda`).

### 6. Compare schemes

    compare_ensembles_jcf()

Prints the comparison table and writes `ensemble_scheme_comparison.csv`
and `jcf_weight_trajectories.csv`.

**Check the output says** `CHECK PASSED`**.** $\lambda = 1$ is
algebraically identical to calibration-only weighting, so it must
reproduce step 3 exactly. If it warns instead, stop and fix that before
reading anything into the other $\lambda$ values.

### 7. JCF plots

    plot_forecasts_jcf()

Figures go to `plots_jcf/` (separate from `plots/`), with per-window
figures in a `plots_jcf/JCF-lam050/` subfolder.

# Part B - Adapting to your own dataset

### 1. Add the data

Put a two-column text file (time, incidence) in `input/`.

### 2. Write one options file per model

Copy `options_forecast_GLM_ebola_scenario1.m` and edit: `cadfilename1`,
`caddisease`, `model.fc` / `model.name`, `params.label`, `params.LB`,
`params.UB`, `params.initial`, `vars.initial` (first observation),
`windowsize1`, `forecastingperiod`.

Parameter bounds matter --- a bound set for Ebola will not suit a
different outbreak scale.

### 3. Write a batch script

Copy `batch_GLM_ebola_scenario1.m`, point it at your data file, options
function, calibration length, horizon and shift.

### 4. Update the config block in each analysis script

The same settings appear at the top of five files. **They must agree.**

  -----------------------------------------------------------------------
  Setting                             Files
  ----------------------------------- -----------------------------------
  `dataset_file`, `calib_period`,     `collect_metrics_simple`,
  `forecast_horizon`, `n_windows`,    `create_ensemble_simple`,
  `window_shift`                      `create_ensemble_jcf`,
                                      `plot_forecasts_simple`,
                                      `plot_forecasts_jcf`

  `model_labels` (short names) and    `create_ensemble_simple`,
  `model_names` (`model.name` from    `create_ensemble_jcf`,
  the options file)                   `plot_forecasts_*`
  -----------------------------------------------------------------------

`model_labels` must match the `Model` column written by
`collect_metrics_simple`; `model_names` must match `model.name` exactly,
including the word "model" (e.g. `'GLM model'`).

### 5. Check the window design

-   $H \leq s$ **is required.** JCF uses the *preceding* window's
    forecast score, which is only fully observed once $H \leq s$.
    `create_ensemble_jcf` errors out otherwise.
-   `n_windows` must not exceed what the series supports:
    $\left\lfloor (N - C - H + 1)/s \right\rfloor$. The batch scripts
    compute this from the data; the analysis scripts do not, so set it
    by hand to match.

### 6. Run Part A steps 1--7 on your data.

# Quick reference

  -------------------------------------------------------------------------------
  Step   Command                               Key output
  ------ ------------------------------------- ----------------------------------
  1      `batch_*`                             `output/` model fits

  2      `collect_metrics_simple()`            `forecast_metrics_summary.mat`

  3      `create_ensemble_simple()`            `ensemble_results.csv`

  4      `plot_model_fit_simple()`,            `plots/`
         `plot_forecasts_simple()`             

  5      `run_jcf_sweep()`                     `ensemble_results_JCF-*.csv`

  6      `compare_ensembles_jcf()`             `ensemble_scheme_comparison.csv`

  7      `plot_forecasts_jcf()`                `plots_jcf/`
  -------------------------------------------------------------------------------

**If you re-run only the weighting**, repeat steps 5--7 only. Steps 1--2
depend on the model fits, not on the weights.

**Two things to remember when reading results.** Window 1 is identical
under JCF and calibration-only weighting by construction. Only windows
2+ tell you whether JCF helped, which is why the comparison table
reports a separate windows-2+ column. And the stacked weight-trajectory
figure is the most informative plot: flat bars mean JCF is doing little,
shifting bars mean it is adapting as the epidemic turns.

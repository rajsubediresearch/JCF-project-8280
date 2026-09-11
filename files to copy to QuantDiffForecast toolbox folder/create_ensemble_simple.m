% =========================================================================
% CREATE ENSEMBLE FORECASTS - SIMPLIFIED VERSION
% For Teaching/Demonstration Purposes
% =========================================================================
% This script builds a single ensemble forecast per rolling window using
% the exact weighted empirical LINEAR POOL:
%
%       getensemble_linearpool_from_curves.m
%
% This is now the only ensemble approach in the pipeline. The earlier
% "mixture of trajectories" route (getensemblemodel.m) has been retired,
% so there is no longer an M1 / M2 distinction anywhere in the outputs.
%
% Weighting:
%   'inverse_wis' (default)  w_i = (1/WIS_i) / sum_j (1/WIS_j)
%                            where WIS_i is the calibration-period WIS of
%                            model i, read from forecast_metrics_summary.mat
%   'equal'                  w_i = 1/K
%
% Two ensembles are formed per window, from the SAME weights:
%   (a) process curves      (forecast_model1)  -> ensemble median
%   (b) predictive curves   (forecast_model12) -> median, 95% PI, WIS
% getensemble_linearpool_from_curves.m documents that it must be called
% separately for the two curve types. Doing so also means the ensemble
% metrics are computed exactly the way Run_Forecasting_ODEModel computes
% them for the individual models (MAE/MSE from the process median, PI
% coverage and WIS from the predictive curves), so the ensemble row is
% directly comparable to the individual model rows.
%
% Requires toolbox functions:
%   - getensemble_linearpool_from_curves.m
%   - computeWIS.m
%
% Run collect_metrics_simple() first (it writes forecast_metrics_summary.mat).
%
% =========================================================================

function create_ensemble_simple()

    fprintf('\n========================================================\n');
    fprintf('   ENSEMBLE FORECAST CREATION (LINEAR POOL)\n');
    fprintf('========================================================\n\n');

    % =====================================================================
    % CONFIGURATION - Modify these for your dataset
    % =====================================================================

    config.dataset_file     = 'ebola-scenario-1';
    config.calib_period     = 20;            % Calibration window (weeks)
    config.forecast_horizon = 4;             % Forecast horizon (weeks)
    config.n_windows        = 5;             % Number of forecast windows
    config.window_shift     = 4;             % Weeks between window starts
    config.output_folder    = './output';    % Where model outputs are stored

    % Models to combine. Add or remove entries here to change the ensemble;
    % nothing below is hard-coded to three models.
    %   model_labels : short name, must match the Model column written by
    %                  collect_metrics_simple()
    %   model_names  : model.name as set in the options_forecast_* file,
    %                  i.e. the string that appears in the output filenames
    config.model_labels = {'GLM', 'Richards', 'Gompertz'};
    config.model_names  = {'GLM model', 'Richards model', 'Gompertz model'};

    config.weighting = 'inverse_wis';   % 'inverse_wis' or 'equal'
    config.fit_index = 1;               % Which fitted variable to ensemble
    config.ens_seed  = 1;               % Seed for the trajectory sample
    config.ens_nens  = [];              % [] = as many draws as each model
                                        % supplies; larger reduces Monte
                                        % Carlo noise in the ensemble WIS
    config.ens_label = 'LinearPool';    % Appears in the output filenames

    % ---- normalise the dataset name (filenames always carry the .txt) ----
    dataset_stem = config.dataset_file;
    if endsWith(lower(dataset_stem), '.txt')
        dataset_stem = dataset_stem(1:end-4);
    end
    config.dataset_stem = dataset_stem;

    K = numel(config.model_labels);
    if numel(config.model_names) ~= K
        error('config.model_labels and config.model_names must be the same length.');
    end

    H = config.forecast_horizon;
    T = config.calib_period + H;

    fprintf('Configuration:\n');
    fprintf('  Dataset: %s\n', config.dataset_stem);
    fprintf('  Models (%d): %s\n', K, strjoin(config.model_labels, ', '));
    fprintf('  Calibration period: %d weeks\n', config.calib_period);
    fprintf('  Forecast horizon: %d weeks\n', H);
    fprintf('  Windows: %d (shift = %d weeks)\n', config.n_windows, config.window_shift);
    fprintf('  Weighting: %s\n', config.weighting);
    fprintf('  Ensemble: exact weighted empirical linear pool\n\n');

    if ~exist(config.output_folder, 'dir')
        error('Output folder not found: %s', config.output_folder);
    end

    % =====================================================================
    % LOAD DATA
    % =====================================================================

    fprintf('Loading data...\n');

    data = load(fullfile('.', 'input', [config.dataset_stem '.txt']));
    fprintf('  Loaded observed data: %d time points\n', size(data, 1));

    if ~exist('forecast_metrics_summary.mat', 'file')
        error('Metrics file not found. Run collect_metrics_simple() first.');
    end
    S = load('forecast_metrics_summary.mat', 'results');
    results = S.results;
    fprintf('  Loaded metrics from %d model-window combinations\n\n', height(results));

    % =====================================================================
    % PREALLOCATE RESULTS
    % =====================================================================

    nW = config.n_windows;

    out.Window    = (1:nW)';
    out.tstart    = zeros(nW, 1);
    out.W         = nan(nW, K);
    out.WIS_calib = nan(nW, 1);
    out.MAE       = nan(nW, 1);
    out.MSE       = nan(nW, 1);
    out.Coverage  = nan(nW, 1);
    out.WIS       = nan(nW, 1);
    out.ok        = false(nW, 1);

    % =====================================================================
    % PROCESS EACH WINDOW
    % =====================================================================

    fprintf('Creating ensemble forecasts...\n\n');

    for w = 1:nW

        tstart = 1 + (w - 1) * config.window_shift;
        out.tstart(w) = tstart;

        fprintf('[%d/%d] Window %d (tstart = %d)\n', w, nW, w, tstart);

        % -----------------------------------------------------------------
        % 1. CHECK THE WINDOW FITS INSIDE THE OBSERVED SERIES
        % -----------------------------------------------------------------

        if tstart + T - 1 > size(data, 1)
            warning('  Window %d extends past the end of the data, skipping.', w);
            continue;
        end

        % -----------------------------------------------------------------
        % 2. CALCULATE WEIGHTS
        % -----------------------------------------------------------------

        [weights, wis_calib, wok] = compute_weights(results, w, config);
        if ~wok
            warning('  Could not build weights for window %d, skipping.', w);
            continue;
        end

        for m = 1:K
            fprintf('    %-12s calibration WIS = %8.3f   weight = %.4f\n', ...
                config.model_labels{m}, wis_calib(m), weights(m));
        end

        out.W(w, :) = weights(:)';

        % -----------------------------------------------------------------
        % 3. LOAD MODEL CURVES (calibration + forecast rows)
        %    process    = forecast_model1  (no observation error)
        %    predictive = forecast_model12 (observation error added)
        % -----------------------------------------------------------------

        curves_proc = cell(1, K);
        curves_pred = cell(1, K);
        missing = false;

        for m = 1:K
            [curves_proc{m}, ok1] = load_curves(config.model_names{m}, tstart, ...
                                                config, 'forecast_model1', T);
            [curves_pred{m}, ok2] = load_curves(config.model_names{m}, tstart, ...
                                                config, 'forecast_model12', T);
            if ~ok1 || ~ok2
                warning('  Missing curves for %s in window %d.', ...
                    config.model_labels{m}, w);
                missing = true;
            end
        end

        if missing
            warning('  Skipping window %d.', w);
            continue;
        end

        curves_proc = harmonize_columns(curves_proc, config.ens_seed, 'process');
        curves_pred = harmonize_columns(curves_pred, config.ens_seed, 'predictive');

        fprintf('    Loaded %d process and %d predictive trajectories per model\n', ...
            size(curves_proc{1}, 2), size(curves_pred{1}, 2));

        % -----------------------------------------------------------------
        % 4. BUILD THE LINEAR-POOL ENSEMBLE
        %    getensemble_linearpool_from_curves(weights, A1, ..., AK, ...)
        %    returns [Qens, forecast1, Qmodels, meta, curvesEns] where
        %      forecast1 : T x 3 [median, LB 2.5%, UB 97.5%]
        %      curvesEns : T x Nens sampled trajectories (needed by
        %                  computeWIS, which expects draws and not quantiles)
        % -----------------------------------------------------------------

        fprintf('    Pooling...\n');

        % (a) process curves -> ensemble median used for MAE / MSE
        [~, fc_proc] = getensemble_linearpool_from_curves( ...
            weights, curves_proc{:});

        % (b) predictive curves -> median, 95% PI, and the draws for WIS
        if isempty(config.ens_nens)
            [~, fc_pred, ~, ~, ens_draws] = getensemble_linearpool_from_curves( ...
                weights, curves_pred{:}, 'Seed', config.ens_seed);
        else
            [~, fc_pred, ~, ~, ens_draws] = getensemble_linearpool_from_curves( ...
                weights, curves_pred{:}, 'Seed', config.ens_seed, ...
                'Nens', config.ens_nens);
        end

        % -----------------------------------------------------------------
        % 5. EVALUATE
        %    datalatest must start at the window start, exactly as
        %    Run_Forecasting_ODEModel passes data(tstart:end,:) to computeWIS.
        % -----------------------------------------------------------------

        datalatest = data(tstart:end, :);
        data_calib = datalatest(1:config.calib_period, :);

        [WISC_ens, WISFS_ens] = computeWIS(data_calib, datalatest, ens_draws, H);

        if isempty(WISFS_ens)
            warning('  computeWIS returned no forecast WIS for window %d, skipping.', w);
            continue;
        end

        observed   = datalatest(config.calib_period + (1:H), 2);
        median_fc  = fc_proc(config.calib_period + (1:H), 1);
        LB         = fc_pred(config.calib_period + (1:H), 2);
        UB         = fc_pred(config.calib_period + (1:H), 3);

        perf = compute_metrics(observed, median_fc, LB, UB, WISFS_ens);

        out.WIS_calib(w) = WISC_ens;
        out.MAE(w)       = perf.MAE;
        out.MSE(w)       = perf.MSE;
        out.Coverage(w)  = perf.Coverage;
        out.WIS(w)       = perf.WIS;
        out.ok(w)        = true;

        fprintf('    Ensemble: MAE=%.2f, MSE=%.2f, Coverage=%.1f%%, WIS=%.2f\n', ...
            perf.MAE, perf.MSE, perf.Coverage, perf.WIS);

        % -----------------------------------------------------------------
        % 6. WRITE OUTPUT FILES FOR THIS WINDOW
        %    Same column layout as the individual-model files, so the
        %    plotting scripts can read them the same way.
        % -----------------------------------------------------------------

        time_vec = datalatest(1:T, 1);
        obs_vec  = datalatest(1:T, 2);

        fc_csv = table(time_vec, obs_vec, fc_pred(:,1), fc_pred(:,2), fc_pred(:,3), ...
            'VariableNames', {'time', 'data', 'median', 'LB', 'UB'});
        writetable(fc_csv, fullfile(config.output_folder, sprintf( ...
            'Forecast-Ensemble-%s-tstart-%d-calibrationperiod-%d-horizon-%d.csv', ...
            config.ens_label, tstart, config.calib_period, H)));

        perf_csv = table((1:H)', perf.MAE_h, perf.MSE_h, perf.Coverage_h, WISFS_ens(:,2), ...
            'VariableNames', {'forecasting_horizon', 'MAE', 'MSE', 'Coverage 95%PI', 'WIS'});
        writetable(perf_csv, fullfile(config.output_folder, sprintf( ...
            'performance-forecasting-Ensemble-%s-tstart-%d-calibrationperiod-%d-horizon-%d.csv', ...
            config.ens_label, tstart, config.calib_period, H)));

        fprintf('    Saved forecast and performance CSVs\n\n');
    end

    % =====================================================================
    % ASSEMBLE RESULTS TABLE
    % =====================================================================

    keep = out.ok;
    if ~any(keep)
        error('No ensemble forecasts were produced. Check the output folder and the configuration.');
    end

    results_ensemble = table(out.Window(keep), out.tstart(keep), ...
        'VariableNames', {'Window', 'tstart'});

    for m = 1:K
        results_ensemble.(weight_varname(config.model_labels{m})) = out.W(keep, m);
    end

    results_ensemble.WIS_calib = out.WIS_calib(keep);
    results_ensemble.MAE       = out.MAE(keep);
    results_ensemble.MSE       = out.MSE(keep);
    results_ensemble.Coverage  = out.Coverage(keep);
    results_ensemble.WIS       = out.WIS(keep);

    % =====================================================================
    % SUMMARY
    % =====================================================================

    fprintf('========================================================\n');
    fprintf('   ENSEMBLE PERFORMANCE SUMMARY\n');
    fprintf('========================================================\n\n');

    fprintf('Linear pool - averaged across %d window(s):\n', sum(keep));
    fprintf('  MAE:      %.2f\n', mean(results_ensemble.MAE));
    fprintf('  MSE:      %.2f\n', mean(results_ensemble.MSE));
    fprintf('  Coverage: %.1f%%\n', mean(results_ensemble.Coverage));
    fprintf('  WIS:      %.2f\n\n', mean(results_ensemble.WIS));

    fprintf('Mean weight by model:\n');
    for m = 1:K
        fprintf('  %-12s %.4f\n', config.model_labels{m}, mean(out.W(keep, m)));
    end
    fprintf('\n');

    % =====================================================================
    % SAVE RESULTS
    % =====================================================================

    save('ensemble_results.mat', 'results_ensemble', 'config');
    fprintf('Saved: ensemble_results.mat\n');

    writetable(results_ensemble, 'ensemble_results.csv');
    fprintf('Saved: ensemble_results.csv\n\n');

    fprintf('========================================================\n');
    fprintf('   DETAILED RESULTS\n');
    fprintf('========================================================\n\n');

    disp(results_ensemble);

    fprintf('Done!\n\n');
end

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function [weights, wis_calib, ok] = compute_weights(results, w, config)
    % Build the model weights for one window.

    K = numel(config.model_labels);
    wis_calib = nan(K, 1);
    weights = nan(K, 1);
    ok = false;

    for m = 1:K
        idx = (results.Window == w) & strcmp(results.Model, config.model_labels{m});
        if ~any(idx)
            warning('  No metrics row for %s in window %d.', config.model_labels{m}, w);
            return;
        end
        v = results.WIS_calib(idx);
        wis_calib(m) = v(1);
    end

    switch lower(config.weighting)

        case 'equal'
            weights = ones(K, 1) / K;

        case 'inverse_wis'
            if any(~isfinite(wis_calib)) || any(wis_calib <= 0)
                warning('  Calibration WIS is missing or non-positive; cannot use inverse-WIS weights.');
                return;
            end
            inv_wis = 1 ./ wis_calib;
            weights = inv_wis / sum(inv_wis);

        otherwise
            error('Unknown config.weighting: %s', config.weighting);
    end

    % getensemble_linearpool_from_curves requires the weights to sum to one
    % within 1e-8, so renormalise defensively after the division above.
    weights = weights / sum(weights);
    ok = true;
end

function [curves, found] = load_curves(model_name, tstart, config, field, T)
    % Load bootstrap trajectories for one model and window.
    %
    % Returns a T x M matrix covering the calibration period followed by the
    % forecast horizon, for the fitted variable given by config.fit_index.

    curves = [];
    found = false;

    % '-tstart-%d-tend-' pins the window exactly. A bare '*tstart-1*' would
    % also match tstart-13, tstart-17, and so on.
    pattern = sprintf(['Forecast-ODEModel-%s.txt-model_name-%s-fit_index-%d-' ...
                       '*-tstart-%d-tend-*-calibrationperiod-%d-*.mat'], ...
        config.dataset_stem, model_name, config.fit_index, tstart, config.calib_period);

    files = dir(fullfile(config.output_folder, pattern));

    if isempty(files)
        return;
    end

    md = load(fullfile(config.output_folder, files(1).name));

    if ~isfield(md, field)
        warning('  %s not present in %s', field, files(1).name);
        return;
    end

    A = md.(field);

    % The saved matrices stack one T-row block per fitted variable, in the
    % order given by vars.fit_index. Select the block we asked for.
    blk = 1;
    if isfield(md, 'vars') && isfield(md.vars, 'fit_index')
        hit = find(md.vars.fit_index == config.fit_index, 1, 'first');
        if ~isempty(hit)
            blk = hit;
        end
    end

    rows = (blk - 1) * T + (1:T);

    if size(A, 1) < rows(end)
        warning('  %s in %s has too few rows (%d) for a %d-row block.', ...
            field, files(1).name, size(A, 1), T);
        return;
    end

    curves = A(rows, :);
    found = true;
end

function C = harmonize_columns(C, seed, tag)
    % All matrices passed to the linear pool must be the same size T x M.
    % If the models were run with different numbers of bootstrap replicates,
    % subsample every model down to the smallest column count.

    n = cellfun(@(x) size(x, 2), C);

    if all(n == n(1))
        return;
    end

    nmin = min(n);
    warning('  Unequal %s draw counts across models (%s); subsampling all to %d.', ...
        tag, strjoin(arrayfun(@(x) num2str(x), n, 'UniformOutput', false), ', '), nmin);

    stream = RandStream('mt19937ar', 'Seed', seed);
    for i = 1:numel(C)
        if size(C{i}, 2) > nmin
            pick = sort(randperm(stream, size(C{i}, 2), nmin));
            C{i} = C{i}(:, pick);
        end
    end
end

function metrics = compute_metrics(observed, median_fc, LB, UB, WISFS)
    % Forecast metrics, aggregated the way Run_Forecasting_ODEModel does it:
    % each horizon h reports the mean over horizons 1..h, and the summary
    % value is the mean across those rows. Matching this keeps the ensemble
    % numbers comparable with the individual model numbers collected by
    % collect_metrics_simple().
    %
    % Inputs:
    %   observed  - observed values over the forecast period (H x 1)
    %   median_fc - ensemble median from the PROCESS curves (H x 1)
    %   LB, UB    - 95% bounds from the PREDICTIVE curves (H x 1)
    %   WISFS     - computeWIS output (H x 2: [horizon, WIS])

    H = numel(observed);

    metrics.MAE_h      = zeros(H, 1);
    metrics.MSE_h      = zeros(H, 1);
    metrics.Coverage_h = zeros(H, 1);

    inside = (observed >= LB) & (observed <= UB);

    for h = 1:H
        err = median_fc(1:h) - observed(1:h);
        metrics.MAE_h(h)      = mean(abs(err));
        metrics.MSE_h(h)      = mean(err .^ 2);
        metrics.Coverage_h(h) = 100 * mean(inside(1:h));
    end

    metrics.MAE      = mean(metrics.MAE_h);
    metrics.MSE      = mean(metrics.MSE_h);
    metrics.Coverage = mean(metrics.Coverage_h);
    metrics.WIS      = mean(WISFS(:, 2));
end

function name = weight_varname(label)
    % Turn a model label into a valid table variable name, e.g. 'w_GLM'.
    name = ['w_' regexprep(label, '\W', '_')];
end

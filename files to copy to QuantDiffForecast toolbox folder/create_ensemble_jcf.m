% =========================================================================
% CREATE ENSEMBLE FORECASTS - JOINT CALIBRATION-FORECASTING (JCF) WEIGHTS
% For Teaching/Demonstration Purposes
% =========================================================================
% Builds one ensemble forecast per rolling window using the exact weighted
% empirical LINEAR POOL:
%
%       getensemble_linearpool_from_curves.m
%
% This is the same ensemble construction used by create_ensemble_simple.m.
% The ONLY difference is how the model weights are derived.
%
% JCF weighting (raw blend):
%
%   S_{i,w} = WIS_calib_{i,1}                                   for w = 1
%   S_{i,w} = lambda * WIS_calib_{i,w}
%             + (1 - lambda) * WIS_forecast_{i,w-1}             for w >= 2
%
%   weight_{i,w} = (1 / S_{i,w}) / sum_j (1 / S_{j,w})
%
% i.e. the current window's calibration WIS blended with the IMMEDIATELY
% PRECEDING window's forecast WIS. Window 1 has no preceding forecast, so it
% falls back to calibration WIS alone - which means window 1 is identical to
% the calibration-only baseline by construction, and only windows 2..W carry
% any signal about whether JCF helps.
%
% lambda = 1 reproduces the calibration-only baseline exactly. Run that first
% as a correctness check against create_ensemble_simple.m before trusting any
% interior lambda.
%
% ALL OUTPUT FILENAMES ARE TAGGED WITH THE SCHEME AND LAMBDA, so this script
% never overwrites the calibration-only pipeline outputs.
%
% Requires toolbox functions:
%   - getensemble_linearpool_from_curves.m
%   - computeWIS.m
%
% Run the batch_* scripts and collect_metrics_simple() first. Neither needs to
% be repeated for JCF - both JCF inputs (WIS_calib and WIS_forecast) are
% already stored in forecast_metrics_summary.mat.
%
% =========================================================================

function create_ensemble_jcf(overrides)

    if nargin < 1
        overrides = struct();
    end

    fprintf('\n========================================================\n');
    fprintf('   ENSEMBLE FORECAST CREATION (JCF WEIGHTS)\n');
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

    % Models to combine. Add or remove entries here; nothing below is
    % hard-coded to three models.
    %   model_labels : must match the Model column from collect_metrics_simple
    %   model_names  : model.name from the options_forecast_* file
    config.model_labels = {'GLM', 'Richards', 'Gompertz'};
    config.model_names  = {'GLM model', 'Richards model', 'Gompertz model'};

    % JCF mixing parameter. Either a scalar in [0,1] or the string 'length'.
    %   1.0      = calibration WIS only (reproduces the baseline pipeline)
    %   0.5      = equal weight to the two periods
    %   0.0      = preceding-window forecast WIS only
    %   'length' = weight each period by its number of time points,
    %              lambda = C / (C + H). Because WIS_calib and WIS_forecast
    %              are already per-point means, this is the same as pooling
    %              the calibration and forecast points into one sample.
    %              It is a useful non-arbitrary reference point, but note
    %              that it treats an in-sample and an out-of-sample point as
    %              equally informative, which is the opposite of JCF's
    %              premise. With C=20, H=4 it gives lambda = 0.833, i.e.
    %              close to calibration-only. It is also dataset-dependent,
    %              so it is not comparable across different window designs.
    config.lambda = 0.5;

    config.fit_index = 1;               % Which fitted variable to ensemble
    config.ens_seed  = 1;               % Seed for the trajectory sample
    config.ens_nens  = [];              % [] = as many draws as each model
                                        % supplies; larger reduces Monte
                                        % Carlo noise in the ensemble WIS

    % ---- apply any caller overrides (used by a lambda sweep driver) ------
    fn = fieldnames(overrides);
    for k = 1:numel(fn)
        config.(fn{k}) = overrides.(fn{k});
    end

    % ---- derived settings ------------------------------------------------
    H = config.forecast_horizon;
    T = config.calib_period + H;

    % Resolve lambda. 'length' means weight each period by its duration.
    if ischar(config.lambda) || (isstring(config.lambda) && isscalar(config.lambda))
        if ~strcmpi(char(config.lambda), 'length')
            error('config.lambda must be a scalar in [0,1] or the string ''length''.');
        end
        config.lambda_mode = 'length';
        config.lambda = config.calib_period / (config.calib_period + H);
    else
        if ~isnumeric(config.lambda) || ~isscalar(config.lambda) || ...
                config.lambda < 0 || config.lambda > 1
            error('config.lambda must be a scalar in [0,1] or the string ''length''.');
        end
        config.lambda_mode = 'fixed';
    end

    % Scheme tag carried into every output filename. The length-weighted run
    % gets its own tag so it cannot collide with a manually chosen lambda
    % that happens to round to the same value.
    if strcmp(config.lambda_mode, 'length')
        config.ens_label = 'JCF-lamlen';
    else
        config.ens_label = sprintf('JCF-lam%03d', round(config.lambda * 100));
    end

    dataset_stem = config.dataset_file;
    if endsWith(lower(dataset_stem), '.txt')
        dataset_stem = dataset_stem(1:end-4);
    end
    config.dataset_stem = dataset_stem;

    K = numel(config.model_labels);
    if numel(config.model_names) ~= K
        error('config.model_labels and config.model_names must be the same length.');
    end

    % =====================================================================
    % NO-LEAKAGE GUARD
    % =====================================================================
    % Window w-1's forecast period ends at t_{w-1} + C + H - 1.
    % Window w's calibration period ends at t_{w-1} + s + C - 1.
    % For the preceding forecast score to be fully observed by the time
    % window w is issued we need H <= s. Otherwise JCF would be scoring
    % against data that is not yet available in real time.

    if H > config.window_shift
        error(['JCF requires forecast_horizon <= window_shift (H=%d, s=%d).\n' ...
               'With H > s the preceding window''s forecast period is not yet ' ...
               'fully observed when the current window is issued, so the ' ...
               'evaluation would be optimistically biased.'], H, config.window_shift);
    end

    fprintf('Configuration:\n');
    fprintf('  Dataset: %s\n', config.dataset_stem);
    fprintf('  Models (%d): %s\n', K, strjoin(config.model_labels, ', '));
    fprintf('  Calibration period: %d weeks\n', config.calib_period);
    fprintf('  Forecast horizon: %d weeks (shift = %d)\n', H, config.window_shift);
    fprintf('  Windows: %d\n', config.n_windows);
    if strcmp(config.lambda_mode, 'length')
        fprintf('  Weighting: JCF raw blend, lambda = C/(C+H) = %d/%d = %.4f\n', ...
            config.calib_period, T, config.lambda);
        fprintf('    (length-weighted: equivalent to pooling the calibration and\n');
        fprintf('     forecast points into a single sample)\n');
    else
        fprintf('  Weighting: JCF raw blend, lambda = %.2f\n', config.lambda);
    end
    if abs(config.lambda - 1) < 1e-12
        fprintf('    (lambda = 1: equivalent to the calibration-only baseline)\n');
    end
    fprintf('  Ensemble: exact weighted empirical linear pool\n');
    fprintf('  Output tag: %s\n\n', config.ens_label);

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

    if ~ismember('WIS_forecast', fieldnames_of(results))
        error(['forecast_metrics_summary.mat has no WIS_forecast column. ' ...
               'JCF needs it; re-run collect_metrics_simple().']);
    end

    % =====================================================================
    % PREALLOCATE RESULTS
    % =====================================================================

    nW = config.n_windows;

    out.Window     = (1:nW)';
    out.tstart     = zeros(nW, 1);
    out.W          = nan(nW, K);   % weights
    out.WIScal_in  = nan(nW, K);   % current-window calibration WIS input
    out.WISfor_in  = nan(nW, K);   % preceding-window forecast WIS input
    out.WIS_calib  = nan(nW, 1);   % ensemble's own calibration WIS
    out.MAE        = nan(nW, 1);
    out.MSE        = nan(nW, 1);
    out.Coverage   = nan(nW, 1);
    out.WIS        = nan(nW, 1);
    out.ok         = false(nW, 1);

    % =====================================================================
    % PROCESS EACH WINDOW
    % =====================================================================

    fprintf('Creating JCF ensemble forecasts...\n\n');

    for w = 1:nW

        tstart = 1 + (w - 1) * config.window_shift;
        out.tstart(w) = tstart;

        fprintf('[%d/%d] Window %d (tstart = %d)\n', w, nW, w, tstart);

        if tstart + T - 1 > size(data, 1)
            warning('  Window %d extends past the end of the data, skipping.', w);
            continue;
        end

        % -----------------------------------------------------------------
        % 1. JCF WEIGHTS
        % -----------------------------------------------------------------

        [weights, wcal, wfor, wok] = compute_jcf_weights(results, w, config);
        if ~wok
            warning('  Could not build JCF weights for window %d, skipping.', w);
            continue;
        end

        if w == 1
            fprintf('    Window 1: no preceding forecast, using calibration WIS only\n');
            for m = 1:K
                fprintf('    %-12s WIScal = %8.3f                      weight = %.4f\n', ...
                    config.model_labels{m}, wcal(m), weights(m));
            end
        else
            for m = 1:K
                fprintf('    %-12s WIScal = %8.3f  WISfor(w-1) = %8.3f  weight = %.4f\n', ...
                    config.model_labels{m}, wcal(m), wfor(m), weights(m));
            end
        end

        out.W(w, :)         = weights(:)';
        out.WIScal_in(w, :) = wcal(:)';
        out.WISfor_in(w, :) = wfor(:)';

        % -----------------------------------------------------------------
        % 2. LOAD MODEL CURVES
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
        % 3. BUILD THE LINEAR-POOL ENSEMBLE
        % -----------------------------------------------------------------

        fprintf('    Pooling...\n');

        [~, fc_proc] = getensemble_linearpool_from_curves( ...
            weights, curves_proc{:});

        if isempty(config.ens_nens)
            [~, fc_pred, ~, ~, ens_draws] = getensemble_linearpool_from_curves( ...
                weights, curves_pred{:}, 'Seed', config.ens_seed);
        else
            [~, fc_pred, ~, ~, ens_draws] = getensemble_linearpool_from_curves( ...
                weights, curves_pred{:}, 'Seed', config.ens_seed, ...
                'Nens', config.ens_nens);
        end

        % -----------------------------------------------------------------
        % 4. EVALUATE
        % -----------------------------------------------------------------

        datalatest = data(tstart:end, :);
        data_calib = datalatest(1:config.calib_period, :);

        [WISC_ens, WISFS_ens] = computeWIS(data_calib, datalatest, ens_draws, H);

        if isempty(WISFS_ens)
            warning('  computeWIS returned no forecast WIS for window %d, skipping.', w);
            continue;
        end

        observed  = datalatest(config.calib_period + (1:H), 2);
        median_fc = fc_proc(config.calib_period + (1:H), 1);
        LB        = fc_pred(config.calib_period + (1:H), 2);
        UB        = fc_pred(config.calib_period + (1:H), 3);

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
        % 5. WRITE OUTPUT FILES (tagged with the scheme + lambda)
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
        error('No JCF ensemble forecasts were produced. Check the output folder and configuration.');
    end

    results_jcf = table(out.Window(keep), out.tstart(keep), ...
        'VariableNames', {'Window', 'tstart'});

    results_jcf.lambda = repmat(config.lambda, sum(keep), 1);
    results_jcf.lambda_mode = repmat({config.lambda_mode}, sum(keep), 1);

    for m = 1:K
        results_jcf.(['w_' safe_name(config.model_labels{m})]) = out.W(keep, m);
    end
    for m = 1:K
        results_jcf.(['WIScal_' safe_name(config.model_labels{m})]) = out.WIScal_in(keep, m);
    end
    for m = 1:K
        results_jcf.(['WISforPrev_' safe_name(config.model_labels{m})]) = out.WISfor_in(keep, m);
    end

    results_jcf.WIS_calib = out.WIS_calib(keep);
    results_jcf.MAE       = out.MAE(keep);
    results_jcf.MSE       = out.MSE(keep);
    results_jcf.Coverage  = out.Coverage(keep);
    results_jcf.WIS       = out.WIS(keep);

    % =====================================================================
    % SUMMARY
    % =====================================================================

    fprintf('========================================================\n');
    fprintf('   JCF ENSEMBLE SUMMARY (%s)\n', lambda_label(config));
    fprintf('========================================================\n\n');

    fprintf('Averaged across %d window(s):\n', sum(keep));
    fprintf('  MAE:      %.2f\n', mean(results_jcf.MAE));
    fprintf('  MSE:      %.2f\n', mean(results_jcf.MSE));
    fprintf('  Coverage: %.1f%%\n', mean(results_jcf.Coverage));
    fprintf('  WIS:      %.2f\n\n', mean(results_jcf.WIS));

    % Windows 2..W are the only ones where JCF differs from the baseline
    later = out.ok & (out.Window >= 2);
    if any(later)
        fprintf('Averaged across windows 2..%d only (the JCF-informative windows):\n', nW);
        fprintf('  MAE:      %.2f\n', mean(out.MAE(later)));
        fprintf('  WIS:      %.2f\n\n', mean(out.WIS(later)));
    end

    fprintf('Mean weight by model:\n');
    for m = 1:K
        fprintf('  %-12s %.4f\n', config.model_labels{m}, mean(out.W(keep, m)));
    end
    fprintf('\n');

    % =====================================================================
    % SAVE RESULTS (tagged filenames - never clobbers the baseline)
    % =====================================================================

    mat_name = sprintf('ensemble_results_%s.mat', config.ens_label);
    csv_name = sprintf('ensemble_results_%s.csv', config.ens_label);

    save(mat_name, 'results_jcf', 'config');
    fprintf('Saved: %s\n', mat_name);

    writetable(results_jcf, csv_name);
    fprintf('Saved: %s\n\n', csv_name);

    fprintf('========================================================\n');
    fprintf('   DETAILED RESULTS\n');
    fprintf('========================================================\n\n');

    disp(results_jcf);

    fprintf('Done!\n\n');
end

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function [weights, wcal, wfor, ok] = compute_jcf_weights(results, w, config)
    % JCF weights for one window: raw blend of the current window's
    % calibration WIS and the preceding window's forecast WIS.

    K = numel(config.model_labels);
    wcal = nan(K, 1);
    wfor = nan(K, 1);
    weights = nan(K, 1);
    ok = false;

    lambda = config.lambda;

    for m = 1:K

        % current window calibration WIS
        idx = (results.Window == w) & strcmp(results.Model, config.model_labels{m});
        if ~any(idx)
            warning('  No metrics row for %s in window %d.', config.model_labels{m}, w);
            return;
        end
        v = results.WIS_calib(idx);
        wcal(m) = v(1);

        % preceding window forecast WIS
        if w >= 2
            idxp = (results.Window == (w-1)) & strcmp(results.Model, config.model_labels{m});
            if ~any(idxp)
                warning('  No metrics row for %s in window %d (needed for JCF).', ...
                    config.model_labels{m}, w - 1);
                return;
            end
            vp = results.WIS_forecast(idxp);
            wfor(m) = vp(1);
        end
    end

    if w == 1
        % No preceding forecast: fall back to calibration WIS alone.
        score = wcal;
    else
        score = lambda * wcal + (1 - lambda) * wfor;
    end

    if any(~isfinite(score)) || any(score <= 0)
        warning('  JCF score is missing or non-positive; cannot invert.');
        return;
    end

    inv_score = 1 ./ score;
    weights = inv_score / sum(inv_score);
    weights = weights / sum(weights);   % defensive renormalisation
    ok = true;
end

function [curves, found] = load_curves(model_name, tstart, config, field, T)
    % Load bootstrap trajectories for one model and window.
    % Returns a T x M matrix: calibration period followed by forecast horizon.

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

    % Saved matrices stack one T-row block per fitted variable, in the order
    % given by vars.fit_index. Select the block we asked for.
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
    % If models were run with different bootstrap counts, subsample down.

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
    % Aggregated the way Run_Forecasting_ODEModel does it: horizon h reports
    % the mean over horizons 1..h, and the summary is the mean across rows.
    % Matching this keeps ensemble numbers comparable with model numbers.

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

function s = safe_name(label)
    s = regexprep(label, '\W', '_');
end

function f = fieldnames_of(t)
    % Column names of a table (or fields of a struct), as a cellstr.
    if isstruct(t)
        f = fieldnames(t);
    else
        f = t.Properties.VariableNames;
    end
end

function s = lambda_label(config)
    % Human-readable description of the mixing parameter.
    if isfield(config, 'lambda_mode') && strcmp(config.lambda_mode, 'length')
        s = sprintf('lambda = %.3f, length-weighted', config.lambda);
    else
        s = sprintf('lambda = %.2f', config.lambda);
    end
end

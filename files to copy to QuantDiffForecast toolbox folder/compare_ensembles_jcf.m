% =========================================================================
% COMPARE ENSEMBLE WEIGHTING SCHEMES
% =========================================================================
% Assembles one comparison table across:
%   - each individual model (from forecast_metrics_summary.mat)
%   - the calibration-only ensemble baseline (ensemble_results.mat)
%   - every JCF run found on disk (ensemble_results_JCF-lam*.mat)
%
% Reports mean MAE, MSE, 95% coverage and WIS, plus a WIS skill score
% relative to the calibration-only ensemble:
%
%   Skill = 100 * (1 - WIS_scheme / WIS_baseline)   [positive = better]
%
% Because window 1 is identical under JCF and calibration-only weighting by
% construction, the table also reports WIS over windows 2..W only. That
% column is the honest comparison - the all-window column is diluted by a
% window in which the two schemes cannot differ.
%
% Writes:
%   ensemble_scheme_comparison.csv   - the comparison table
%   jcf_weight_trajectories.csv      - weights by window and lambda
%
% This script only reads; it never overwrites pipeline outputs.
%
% Run after create_ensemble_simple() and run_jcf_sweep().
% =========================================================================

function compare_ensembles_jcf()

    fprintf('\n========================================================\n');
    fprintf('   ENSEMBLE SCHEME COMPARISON\n');
    fprintf('========================================================\n\n');

    config.baseline_file = 'ensemble_results.mat';
    config.jcf_pattern   = 'ensemble_results_JCF-lam*.mat';

    % =====================================================================
    % LOAD
    % =====================================================================

    if ~exist('forecast_metrics_summary.mat', 'file')
        error('forecast_metrics_summary.mat not found. Run collect_metrics_simple() first.');
    end
    S = load('forecast_metrics_summary.mat', 'results');
    results = S.results;

    have_baseline = (exist(config.baseline_file, 'file') == 2);
    if have_baseline
        B = load(config.baseline_file, 'results_ensemble');
        baseline = B.results_ensemble;
        fprintf('Baseline (calibration-only): %s\n', config.baseline_file);
    else
        baseline = [];
        warning('Baseline %s not found; skill scores will be omitted.', config.baseline_file);
    end

    jcf_files = dir(config.jcf_pattern);
    if isempty(jcf_files)
        warning('No JCF result files matching %s. Run run_jcf_sweep() first.', config.jcf_pattern);
    else
        fprintf('Found %d JCF run(s)\n', numel(jcf_files));
    end
    fprintf('\n');

    % =====================================================================
    % ROWS: INDIVIDUAL MODELS
    % =====================================================================

    model_labels = unique(results.Model, 'stable');
    nM = numel(model_labels);

    scheme  = {};
    MAE = []; MSE = []; COV = []; WIS = []; WIS2 = [];

    for m = 1:nM
        idx = strcmp(results.Model, model_labels{m});
        scheme{end+1,1} = model_labels{m}; %#ok<AGROW>
        MAE(end+1,1) = mean(results.MAE_forecast(idx), 'omitnan');       %#ok<AGROW>
        MSE(end+1,1) = mean(results.MSE_forecast(idx), 'omitnan');       %#ok<AGROW>
        COV(end+1,1) = mean(results.Coverage_forecast(idx), 'omitnan');  %#ok<AGROW>
        WIS(end+1,1) = mean(results.WIS_forecast(idx), 'omitnan');       %#ok<AGROW>

        later = idx & (results.Window >= 2);
        WIS2(end+1,1) = mean(results.WIS_forecast(later), 'omitnan');    %#ok<AGROW>
    end

    % =====================================================================
    % ROW: CALIBRATION-ONLY ENSEMBLE BASELINE
    % =====================================================================

    baseline_WIS = NaN;

    if have_baseline
        scheme{end+1,1} = 'Ensemble: calibration-only';
        MAE(end+1,1) = mean(baseline.MAE, 'omitnan');
        MSE(end+1,1) = mean(baseline.MSE, 'omitnan');
        COV(end+1,1) = mean(baseline.Coverage, 'omitnan');
        WIS(end+1,1) = mean(baseline.WIS, 'omitnan');

        lb = baseline.Window >= 2;
        WIS2(end+1,1) = mean(baseline.WIS(lb), 'omitnan');

        baseline_WIS = WIS(end);
    end

    % =====================================================================
    % ROWS: JCF RUNS (sorted by lambda)
    % =====================================================================

    jcf_row_start = numel(scheme) + 1;   % first JCF row in the table

    lam_list = [];
    lam_desc = {};
    jcf_tables = {};

    for f = 1:numel(jcf_files)
        J = load(jcf_files(f).name, 'results_jcf', 'config');
        lam_list(end+1,1) = J.config.lambda;             %#ok<AGROW>
        lam_desc{end+1,1} = lambda_desc(J.config);       %#ok<AGROW>
        jcf_tables{end+1} = J.results_jcf;               %#ok<AGROW>
    end

    [lam_list, order] = sort(lam_list);
    lam_desc = lam_desc(order);
    jcf_tables = jcf_tables(order);

    for f = 1:numel(jcf_tables)
        rj = jcf_tables{f};
        scheme{end+1,1} = sprintf('Ensemble: JCF %s', lam_desc{f});
        MAE(end+1,1) = mean(rj.MAE, 'omitnan');
        MSE(end+1,1) = mean(rj.MSE, 'omitnan');
        COV(end+1,1) = mean(rj.Coverage, 'omitnan');
        WIS(end+1,1) = mean(rj.WIS, 'omitnan');

        lb = rj.Window >= 2;
        WIS2(end+1,1) = mean(rj.WIS(lb), 'omitnan');
    end

    % =====================================================================
    % SKILL SCORES
    % =====================================================================

    if isfinite(baseline_WIS) && baseline_WIS > 0
        Skill = 100 * (1 - WIS / baseline_WIS);
    else
        Skill = nan(size(WIS));
    end

    comparison = table(scheme, MAE, MSE, COV, WIS, WIS2, Skill, ...
        'VariableNames', {'Scheme', 'MAE', 'MSE', 'Coverage', ...
                          'WIS', 'WIS_windows2plus', 'WIS_skill_pct'});

    % =====================================================================
    % PRINT
    % =====================================================================

    fprintf('%-30s | %8s | %8s | %8s | %8s | %8s\n', ...
        'Scheme', 'MAE', 'MSE', 'Cover%', 'WIS', 'Skill%');
    fprintf('%s\n', repmat('-', 1, 88));
    for r = 1:numel(scheme)
        fprintf('%-30s | %8.2f | %8.2f | %8.1f | %8.2f | %8.1f\n', ...
            scheme{r}, MAE(r), MSE(r), COV(r), WIS(r), Skill(r));
    end
    fprintf('\n');
    fprintf('Skill is relative to the calibration-only ensemble (positive = better).\n');
    fprintf('WIS_windows2plus in the CSV excludes window 1, where JCF and the\n');
    fprintf('calibration-only baseline are identical by construction.\n\n');

    if ~isempty(lam_list)
        fprintf('Lambda sensitivity (mean forecast WIS):\n');
        for f = 1:numel(lam_list)
            k = jcf_row_start + f - 1;
            fprintf('  %-22s WIS = %7.3f   (windows 2+: %7.3f)\n', ...
                lam_desc{f}, WIS(k), WIS2(k));
        end
        fprintf('\n');

        % Correctness check: lambda = 1 must match the baseline
        i1 = find(abs(lam_list - 1) < 1e-12, 1);
        if ~isempty(i1) && have_baseline
            k = jcf_row_start + i1 - 1;
            d = abs(WIS(k) - baseline_WIS);
            if d < 1e-8
                fprintf('CHECK PASSED: lambda = 1 reproduces the calibration-only baseline.\n\n');
            else
                warning(['lambda = 1 does NOT match the calibration-only baseline ' ...
                         '(WIS differs by %.3g). Investigate before trusting other lambdas.'], d);
            end
        end
    end

    % =====================================================================
    % WEIGHT TRAJECTORIES
    % =====================================================================

    if ~isempty(jcf_tables)
        traj_lambda = [];
        traj_window = [];
        traj_model  = {};
        traj_weight = [];

        for f = 1:numel(jcf_tables)
            rj = jcf_tables{f};
            vn = varnames_of(rj);
            wcols = vn(startsWith(vn, 'w_'));

            for c = 1:numel(wcols)
                lbl = wcols{c}(3:end);
                for r = 1:height(rj)
                    traj_lambda(end+1,1) = lam_list(f);   %#ok<AGROW>
                    traj_window(end+1,1) = rj.Window(r);  %#ok<AGROW>
                    traj_model{end+1,1}  = lbl;           %#ok<AGROW>
                    traj_weight(end+1,1) = rj.(wcols{c})(r); %#ok<AGROW>
                end
            end
        end

        traj = table(traj_lambda, traj_window, traj_model, traj_weight, ...
            'VariableNames', {'lambda', 'Window', 'Model', 'weight'});

        writetable(traj, 'jcf_weight_trajectories.csv');
        fprintf('Saved: jcf_weight_trajectories.csv\n');
    end

    % =====================================================================
    % SAVE
    % =====================================================================

    writetable(comparison, 'ensemble_scheme_comparison.csv');
    fprintf('Saved: ensemble_scheme_comparison.csv\n\n');

    fprintf('Done!\n\n');
end

function s = lambda_desc(cfg)
    % Short description of a run's mixing parameter.
    if isfield(cfg, 'lambda_mode') && strcmp(cfg.lambda_mode, 'length')
        s = sprintf('lambda=%.3f (length)', cfg.lambda);
    else
        s = sprintf('lambda=%.2f', cfg.lambda);
    end
end

function f = varnames_of(t)
    % Column names of a table (or fields of a struct), as a cellstr.
    if isstruct(t)
        f = fieldnames(t);
    else
        f = t.Properties.VariableNames;
    end
end

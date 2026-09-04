% =========================================================================
% PLOT JCF ENSEMBLE RESULTS
% =========================================================================
% Produces the figure set for the JCF weighting scheme. Everything is
% written to a SEPARATE folder (./plots_jcf by default), so the figures from
% the calibration-only pipeline in ./plots are never overwritten.
%
% Figures produced:
%
%   ./plots_jcf/<tag>/                    (one folder per lambda run found)
%       window<w>_jcf_ensemble            forecast, 95% PI, observed data
%       window<w>_jcf_vs_models           ensemble vs individual model medians
%
%   ./plots_jcf/
%       jcf_weight_trajectories           stacked weights by window, per lambda
%       jcf_lambda_sensitivity            mean forecast WIS vs lambda
%       jcf_scheme_comparison             models vs baseline vs JCF, 4 metrics
%       jcf_metrics_by_window             per-window WIS for every scheme
%
% Per-window figures are produced only for the "focus" lambda by default
% (config.lambda_focus); set config.all_lambdas_per_window = true to produce
% them for every lambda, which is a lot of files.
%
% This script only reads results; it never recomputes or overwrites them.
% Run after create_ensemble_simple(), run_jcf_sweep() and (optionally)
% compare_ensembles_jcf().
% =========================================================================

function plot_forecasts_jcf()

    fprintf('\n========================================================\n');
    fprintf('   PLOTTING JCF ENSEMBLE RESULTS\n');
    fprintf('========================================================\n\n');

    % =====================================================================
    % CONFIGURATION - keep in sync with create_ensemble_jcf()
    % =====================================================================

    config.dataset_file     = 'ebola-scenario-1';
    config.calib_period     = 20;
    config.forecast_horizon = 4;
    config.n_windows        = 5;
    config.window_shift     = 4;
    config.output_folder    = './output';
    config.plots_folder     = './plots_jcf';   % SEPARATE from ./plots

    config.model_labels = {'GLM', 'Richards', 'Gompertz'};

    config.lambda_focus = 0.5;   % lambda used for the per-window figures
    config.all_lambdas_per_window = false;

    config.baseline_file = 'ensemble_results.mat';
    config.jcf_pattern   = 'ensemble_results_JCF-lam*.mat';

    % ---------------------------------------------------------------------

    models = config.model_labels;
    K = numel(models);

    model_colors = num2cell(lines(K), 2);
    ens_color    = [0.85, 0.33, 0.10];   % JCF ensemble (orange)
    base_color   = [0.30, 0.30, 0.30];   % calibration-only baseline (grey)

    if ~exist(config.plots_folder, 'dir')
        mkdir(config.plots_folder);
        fprintf('Created folder: %s\n', config.plots_folder);
    end

    % =====================================================================
    % LOAD RESULTS
    % =====================================================================

    jcf_files = dir(config.jcf_pattern);
    if isempty(jcf_files)
        error(['No JCF result files matching %s. ' ...
               'Run run_jcf_sweep() or create_ensemble_jcf() first.'], config.jcf_pattern);
    end

    lam_list = [];
    lam_short = {};
    jcf_tables = {};
    jcf_tags = {};

    for f = 1:numel(jcf_files)
        J = load(jcf_files(f).name, 'results_jcf', 'config');
        lam_list(end+1,1)  = J.config.lambda;           %#ok<AGROW>
        lam_short{end+1,1} = lambda_short(J.config);    %#ok<AGROW>
        jcf_tables{end+1}  = J.results_jcf;             %#ok<AGROW>
        jcf_tags{end+1}    = J.config.ens_label;        %#ok<AGROW>
    end

    [lam_list, order] = sort(lam_list);
    lam_short  = lam_short(order);
    jcf_tables = jcf_tables(order);
    jcf_tags   = jcf_tags(order);
    nL = numel(lam_list);

    fprintf('Found %d JCF run(s): lambda = %s\n', nL, mat2str(lam_list', 3));

    have_baseline = (exist(config.baseline_file, 'file') == 2);
    if have_baseline
        B = load(config.baseline_file, 'results_ensemble');
        baseline = B.results_ensemble;
        fprintf('Baseline loaded: %s\n', config.baseline_file);
    else
        baseline = [];
        fprintf('Baseline %s not found; baseline curves will be omitted.\n', config.baseline_file);
    end

    have_metrics = (exist('forecast_metrics_summary.mat', 'file') == 2);
    if have_metrics
        S = load('forecast_metrics_summary.mat', 'results');
        results = S.results;
    else
        results = [];
        fprintf('forecast_metrics_summary.mat not found; model rows will be omitted.\n');
    end

    dataset_stem = config.dataset_file;
    if endsWith(lower(dataset_stem), '.txt')
        dataset_stem = dataset_stem(1:end-4);
    end
    config.dataset_stem = dataset_stem;

    % Which lambda runs get per-window figures
    if config.all_lambdas_per_window
        focus_idx = 1:nL;
    else
        [~, focus_idx] = min(abs(lam_list - config.lambda_focus));
        fprintf('Per-window figures for %s\n', jcf_tags{focus_idx});
    end
    fprintf('\n');

    % =====================================================================
    % 1. PER-WINDOW FIGURES
    % =====================================================================

    for fi = focus_idx

        tag = jcf_tags{fi};
        sub = fullfile(config.plots_folder, tag);
        if ~exist(sub, 'dir')
            mkdir(sub);
        end

        fprintf('Per-window figures for %s -> %s\n', tag, sub);

        for w = 1:config.n_windows

            tstart = 1 + (w - 1) * config.window_shift;

            ens_csv = ensemble_csv_path(tag, tstart, config);
            if ~exist(ens_csv, 'file')
                fprintf('  Window %d: ensemble CSV not found, skipping.\n', w);
                continue;
            end

            fc_ens = readtable(ens_csv, 'VariableNamingRule', 'preserve');

            % -------- (a) ensemble forecast on its own -------------------

            fig = figure('Name', sprintf('%s window %d', tag, w), ...
                         'Position', [100, 100, 700, 450], 'Visible', 'off');

            plot_forecast_panel(fc_ens, config, ens_color);
            title(sprintf('JCF Ensemble (%s) - Window %d (tstart = %d)', ...
                lam_short{fi}, w, tstart));

            save_figure(fig, fullfile(sub, sprintf('window%d_jcf_ensemble', w)));

            % -------- (b) ensemble vs individual models ------------------

            fig = figure('Name', sprintf('%s window %d vs models', tag, w), ...
                         'Position', [100, 100, 1200, 500], 'Visible', 'off');

            hold on;

            handles = gobjects(0);
            labels  = {};

            h = fill([fc_ens.time; flipud(fc_ens.time)], ...
                     [fc_ens.LB; flipud(fc_ens.UB)], ...
                     ens_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            handles(end+1) = h;  labels{end+1} = 'JCF 95% PI';

            for m = 1:K
                mcsv = find_forecast_csv(models{m}, tstart, config);
                if isempty(mcsv)
                    continue;
                end
                fc = readtable(mcsv, 'VariableNamingRule', 'preserve');
                h = plot(fc.time, fc.median, '--', 'Color', model_colors{m}, 'LineWidth', 1.5);
                handles(end+1) = h;  labels{end+1} = models{m};
            end

            h = plot(fc_ens.time, fc_ens.median, '-', 'Color', ens_color, 'LineWidth', 2.5);
            handles(end+1) = h;  labels{end+1} = sprintf('JCF (%s)', lam_short{fi});

            h = plot(fc_ens.time, fc_ens.data, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
            handles(end+1) = h;  labels{end+1} = 'Observed';

            xline(fc_ens.time(config.calib_period), 'k--', 'LineWidth', 2);

            xlabel('Time (weeks)');
            ylabel('Cases');
            title(sprintf('JCF Ensemble vs Individual Models - Window %d (tstart = %d)', w, tstart));
            legend(handles, labels, 'Location', 'best');
            set(gca, 'FontSize', 12);
            grid on;  box on;

            save_figure(fig, fullfile(sub, sprintf('window%d_jcf_vs_models', w)));
        end
        fprintf('\n');
    end

    % =====================================================================
    % 2. WEIGHT TRAJECTORIES (stacked bars, one panel per lambda)
    % =====================================================================

    fprintf('Plotting weight trajectories...\n');

    fig = figure('Name', 'JCF weight trajectories', ...
                 'Position', [100, 100, 350 * min(nL, 3), 300 * ceil(nL / 3)], ...
                 'Visible', 'off');

    for f = 1:nL
        subplot(ceil(nL / 3), min(nL, 3), f);

        rj = jcf_tables{f};
        Wm = weights_matrix(rj, models);

        b = bar(rj.Window, Wm, 'stacked');
        for m = 1:K
            b(m).FaceColor = model_colors{m};
        end

        xlabel('Window');
        ylabel('Weight');
        ylim([0, 1]);
        title(lam_short{f});
        if f == 1
            legend(models, 'Location', 'southoutside', 'Orientation', 'horizontal');
        end
        set(gca, 'FontSize', 9);
        box on;
    end

    sgtitle('JCF model weights by window');
    save_figure(fig, fullfile(config.plots_folder, 'jcf_weight_trajectories'));

    % =====================================================================
    % 3. LAMBDA SENSITIVITY
    % =====================================================================

    fprintf('Plotting lambda sensitivity...\n');

    wis_all  = nan(nL, 1);
    wis_late = nan(nL, 1);

    for f = 1:nL
        rj = jcf_tables{f};
        wis_all(f) = mean(rj.WIS, 'omitnan');
        late = rj.Window >= 2;
        wis_late(f) = mean(rj.WIS(late), 'omitnan');
    end

    fig = figure('Name', 'JCF lambda sensitivity', ...
                 'Position', [100, 100, 700, 450], 'Visible', 'off');
    hold on;

    plot(lam_list, wis_all, 'o-', 'Color', ens_color, 'LineWidth', 2, ...
         'MarkerSize', 8, 'MarkerFaceColor', ens_color);
    plot(lam_list, wis_late, 's--', 'Color', ens_color, 'LineWidth', 2, ...
         'MarkerSize', 8);

    if have_baseline
        yline(mean(baseline.WIS, 'omitnan'), '-', 'Calibration-only baseline', ...
              'Color', base_color, 'LineWidth', 1.5);
    end

    xlabel('\lambda   (1 = calibration only,  0 = preceding forecast only)');
    ylabel('Mean forecast WIS');
    title('JCF sensitivity to the mixing parameter (lower is better)');
    legend({'All windows', 'Windows 2+ only'}, 'Location', 'best');
    set(gca, 'FontSize', 11);
    grid on;  box on;

    save_figure(fig, fullfile(config.plots_folder, 'jcf_lambda_sensitivity'));

    % =====================================================================
    % 4. SCHEME COMPARISON (4 metrics)
    % =====================================================================

    if have_metrics

        fprintf('Plotting scheme comparison...\n');

        metric_names = {'MAE', 'MSE', 'Coverage', 'WIS'};
        model_cols   = {'MAE_forecast', 'MSE_forecast', 'Coverage_forecast', 'WIS_forecast'};

        bar_labels = models;
        bar_colors = model_colors;

        if have_baseline
            bar_labels{end+1} = 'Cal-only';
            bar_colors{end+1} = base_color;
        end
        for f = 1:nL
            bar_labels{end+1} = lam_short{f}; %#ok<AGROW>
            bar_colors{end+1} = ens_color;                             %#ok<AGROW>
        end

        nBar = numel(bar_labels);

        fig = figure('Name', 'Scheme comparison', ...
                     'Position', [100, 100, 1100, 420], 'Visible', 'off');

        for i = 1:4
            subplot(1, 4, i);

            vals = nan(nBar, 1);
            k = 0;

            for m = 1:K
                k = k + 1;
                idx = strcmp(results.Model, models{m});
                vals(k) = mean(results.(model_cols{i})(idx), 'omitnan');
            end

            if have_baseline
                k = k + 1;
                vals(k) = mean(baseline.(metric_names{i}), 'omitnan');
            end

            for f = 1:nL
                k = k + 1;
                vals(k) = mean(jcf_tables{f}.(metric_names{i}), 'omitnan');
            end

            b = bar(vals);
            b.FaceColor = 'flat';
            for j = 1:nBar
                b.CData(j, :) = bar_colors{j};
            end

            set(gca, 'XTick', 1:nBar, 'XTickLabel', bar_labels);
            xtickangle(45);
            ylabel(metric_names{i});
            title(metric_names{i});
            set(gca, 'FontSize', 9);
            grid on;
        end

        sgtitle('Individual models vs calibration-only vs JCF ensembles');
        save_figure(fig, fullfile(config.plots_folder, 'jcf_scheme_comparison'));

        % =================================================================
        % 5. WIS BY WINDOW
        % =================================================================

        fprintf('Plotting WIS by window...\n');

        fig = figure('Name', 'WIS by window', ...
                     'Position', [100, 100, 800, 500], 'Visible', 'off');
        hold on;

        leg = {};
        all_markers = {'o-', 's-', 'd-', '^-', 'v-', 'p-', 'h-', '*-'};

        for m = 1:K
            idx = strcmp(results.Model, models{m});
            mk = all_markers{mod(m - 1, numel(all_markers)) + 1};
            plot(results.Window(idx), results.WIS_forecast(idx), mk, ...
                 'Color', model_colors{m}, 'LineWidth', 1.5, 'MarkerSize', 7, ...
                 'MarkerFaceColor', model_colors{m});
            leg{end+1} = models{m}; %#ok<AGROW>
        end

        if have_baseline
            plot(baseline.Window, baseline.WIS, '-', 'Color', base_color, ...
                 'LineWidth', 2.5, 'Marker', 'o', 'MarkerSize', 8, ...
                 'MarkerFaceColor', base_color);
            leg{end+1} = 'Cal-only ensemble';
        end

        for fi2 = focus_idx
            rj = jcf_tables{fi2};
            plot(rj.Window, rj.WIS, '-', 'Color', ens_color, 'LineWidth', 3, ...
                 'Marker', 'o', 'MarkerSize', 9, 'MarkerFaceColor', ens_color);
            leg{end+1} = sprintf('JCF %s', lam_short{fi2}); %#ok<AGROW>
        end

        % Window 1 is identical under both schemes by construction
        xline(1.5, ':', 'JCF differs from here on', 'Color', [0.4 0.4 0.4], ...
              'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');

        xlabel('Window');
        ylabel('Mean forecast WIS');
        title('Forecast WIS by window (lower is better)');
        legend(leg, 'Location', 'best');
        set(gca, 'FontSize', 11);
        grid on;  box on;

        save_figure(fig, fullfile(config.plots_folder, 'jcf_metrics_by_window'));
    end

    fprintf('\n========================================================\n');
    fprintf('   PLOTTING COMPLETE!\n');
    fprintf('========================================================\n');
    fprintf('   All JCF plots saved to: %s/\n', config.plots_folder);
    fprintf('========================================================\n\n');
end

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function Wm = weights_matrix(rj, models)
    % Pull the w_<model> columns out of a JCF results table, in model order.

    n = numel(rj.Window);
    K = numel(models);
    Wm = nan(n, K);

    for m = 1:K
        col = ['w_' regexprep(models{m}, '\W', '_')];
        vn = varnames_of(rj);
        if any(strcmp(vn, col))
            Wm(:, m) = rj.(col);
        end
    end
end

function plot_forecast_panel(fc, config, color)
    % Shared panel: shaded 95% PI, median, observed points, calibration line.

    t_all = fc.time;

    hold on;

    fill([t_all; flipud(t_all)], [fc.LB; flipud(fc.UB)], ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    plot(t_all, fc.median, 'Color', color, 'LineWidth', 2);

    plot(t_all, fc.data, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');

    xline(t_all(config.calib_period), 'k--', 'LineWidth', 1.5);

    xlabel('Time (weeks)');
    ylabel('Cases');
    legend('95% PI', 'Median', 'Observed', 'Location', 'best');

    set(gca, 'FontSize', 10);
    grid on;  box on;
end

function p = ensemble_csv_path(tag, tstart, config)
    p = fullfile(config.output_folder, sprintf( ...
        'Forecast-Ensemble-%s-tstart-%d-calibrationperiod-%d-horizon-%d.csv', ...
        tag, tstart, config.calib_period, config.forecast_horizon));
end

function csv_file = find_forecast_csv(model_name, tstart, config)
    % '-tstart-%d-tend-' pins the window exactly; a bare '*tstart-1*' would
    % also match tstart-13, tstart-17, and so on.

    pattern = sprintf('Forecast-model_name-%s*-tstart-%d-tend-*calibrationperiod-%d-*.csv', ...
        model_name, tstart, config.calib_period);

    files = dir(fullfile(config.output_folder, pattern));

    if isempty(files)
        csv_file = '';
    else
        csv_file = fullfile(config.output_folder, files(1).name);
    end
end

function s = lambda_short(cfg)
    % Short axis/legend label for a run's mixing parameter.
    if isfield(cfg, 'lambda_mode') && strcmp(cfg.lambda_mode, 'length')
        s = sprintf('\\lambda=%.2f (len)', cfg.lambda);
    else
        s = sprintf('\\lambda=%.2f', cfg.lambda);
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

function save_figure(fig, filepath)
    saveas(fig, [filepath '.png']);
    saveas(fig, [filepath '.fig']);
    close(fig);

    fprintf('    Saved: %s.png and %s.fig\n', filepath, filepath);
end

% =========================================================================
% PLOT FORECASTS - SIMPLIFIED VERSION
% For Teaching/Demonstration Purposes
% =========================================================================
% This script creates plots for:
%   1. Individual model forecasts per window
%   2. Linear pool ensemble forecast per window
%   3. Ensemble vs individual models per window
%   4. Performance metrics comparison (models vs ensemble)
%   5. Forecast metrics by window
%
% The pipeline now uses a single ensemble approach (the exact weighted
% empirical linear pool from getensemble_linearpool_from_curves.m), so
% there is no longer an M1 / M2 pair to plot side by side.
%
% Uses existing output files - no recalculation needed.
% Run collect_metrics_simple() and create_ensemble_simple() first.
%
% =========================================================================

function plot_forecasts_simple()

    fprintf('\n========================================================\n');
    fprintf('   PLOTTING FORECASTS\n');
    fprintf('========================================================\n\n');

    % =====================================================================
    % CONFIGURATION - keep in sync with create_ensemble_simple()
    % =====================================================================

    config.dataset_file     = 'ebola-scenario-1';
    config.dataset_name     = 'Ebola Scenario 1';
    config.calib_period     = 20;
    config.forecast_horizon = 4;
    config.n_windows        = 5;
    config.window_shift     = 4;
    config.output_folder    = './output';
    config.plots_folder     = './plots';
    config.ens_label        = 'LinearPool';

    config.model_labels = {'GLM', 'Richards', 'Gompertz'};

    if ~exist(config.plots_folder, 'dir')
        mkdir(config.plots_folder);
        fprintf('Created folder: %s\n\n', config.plots_folder);
    end

    models = config.model_labels;
    K = numel(models);
    H = config.forecast_horizon;

    % One colour per model, generated so this still works with more models
    model_colors = num2cell(lines(K), 2);
    ens_color = [0.85, 0.33, 0.10];   % orange for the ensemble

    % =====================================================================
    % LOAD DATA
    % =====================================================================

    fprintf('Loading data...\n');
    dataset_stem = config.dataset_file;
    if endsWith(lower(dataset_stem), '.txt')
        dataset_stem = dataset_stem(1:end-4);
    end
    data = load(fullfile('.', 'input', [dataset_stem '.txt']));
    fprintf('  Loaded %d time points\n\n', size(data, 1));

    % =====================================================================
    % 1. PLOT INDIVIDUAL MODEL FORECASTS (per window)
    % =====================================================================

    fprintf('Plotting individual model forecasts...\n\n');

    for w = 1:config.n_windows

        tstart = 1 + (w - 1) * config.window_shift;

        fprintf('  Window %d (tstart = %d)\n', w, tstart);

        fig = figure('Name', sprintf('Window %d Forecasts', w), ...
                     'Position', [100, 100, 400 * K, 400], 'Visible', 'off');

        for m = 1:K

            subplot(1, K, m);

            forecast_csv = find_forecast_csv(models{m}, tstart, config);

            if isempty(forecast_csv)
                title(sprintf('%s - Not Found', models{m}));
                continue;
            end

            fc = readtable(forecast_csv, 'VariableNamingRule', 'preserve');

            plot_forecast_panel(fc, config, model_colors{m});
            title(sprintf('%s Model', models{m}));

        end

        sgtitle(sprintf('Window %d (tstart = %d, calib = %d weeks)', ...
            w, tstart, config.calib_period));

        save_figure(fig, fullfile(config.plots_folder, sprintf('window%d_models', w)));

    end

    % =====================================================================
    % 2. PLOT THE ENSEMBLE FORECAST (per window)
    % =====================================================================

    fprintf('\nPlotting ensemble forecasts...\n\n');

    ensemble_files = dir(fullfile(config.output_folder, ...
        sprintf('Forecast-Ensemble-%s-*.csv', config.ens_label)));

    if isempty(ensemble_files)
        fprintf('  Ensemble forecast CSVs not found. Run create_ensemble_simple() first.\n');
    else

        for w = 1:config.n_windows

            tstart = 1 + (w - 1) * config.window_shift;

            ens_csv = ensemble_csv_path(tstart, config);
            if ~exist(ens_csv, 'file')
                fprintf('  Window %d ensemble CSV not found, skipping.\n', w);
                continue;
            end

            fig = figure('Name', sprintf('Window %d Ensemble Forecast', w), ...
                         'Position', [100, 100, 700, 450], 'Visible', 'off');

            fc = readtable(ens_csv, 'VariableNamingRule', 'preserve');

            plot_forecast_panel(fc, config, ens_color);
            title(sprintf('Linear Pool Ensemble - Window %d (tstart = %d)', w, tstart));

            save_figure(fig, fullfile(config.plots_folder, sprintf('window%d_ensemble', w)));
        end

        % =================================================================
        % 3. PLOT ENSEMBLE VS INDIVIDUAL MODELS (per window)
        % =================================================================

        fprintf('\nPlotting ensemble vs models comparison (all windows)...\n\n');

        for w = 1:config.n_windows

            tstart = 1 + (w - 1) * config.window_shift;

            ens_csv = ensemble_csv_path(tstart, config);
            if ~exist(ens_csv, 'file')
                continue;
            end

            fig = figure('Name', sprintf('Window %d - Ensemble vs Models', w), ...
                         'Position', [100, 100, 1200, 500], 'Visible', 'off');

            hold on;

            handles = gobjects(0);
            labels  = {};

            % Ensemble 95% PI first, so the model lines draw on top of it
            fc_ens = readtable(ens_csv, 'VariableNamingRule', 'preserve');

            h = fill([fc_ens.time; flipud(fc_ens.time)], ...
                     [fc_ens.LB; flipud(fc_ens.UB)], ...
                     ens_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            handles(end+1) = h;
            labels{end+1}  = 'Ensemble 95% PI';

            % Individual model medians (dashed)
            for m = 1:K
                forecast_csv = find_forecast_csv(models{m}, tstart, config);
                if isempty(forecast_csv)
                    continue;
                end
                fc = readtable(forecast_csv, 'VariableNamingRule', 'preserve');
                h = plot(fc.time, fc.median, '--', 'Color', model_colors{m}, 'LineWidth', 1.5);
                handles(end+1) = h;
                labels{end+1}  = models{m};
            end

            % Ensemble median
            h = plot(fc_ens.time, fc_ens.median, '-', 'Color', ens_color, 'LineWidth', 2.5);
            handles(end+1) = h;
            labels{end+1}  = 'Ensemble (Linear Pool)';

            % Observed data
            h = plot(fc_ens.time, fc_ens.data, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
            handles(end+1) = h;
            labels{end+1}  = 'Observed';

            xline(fc_ens.time(config.calib_period), 'k--', 'LineWidth', 2);

            xlabel('Time (weeks)');
            ylabel('Cases');
            title(sprintf('Ensemble vs Individual Models - Window %d (tstart = %d)', w, tstart));
            legend(handles, labels, 'Location', 'best');

            set(gca, 'FontSize', 12);
            grid on;
            box on;

            save_figure(fig, fullfile(config.plots_folder, ...
                sprintf('window%d_ensemble_vs_models', w)));

        end

    end

    % =====================================================================
    % 4. PLOT PERFORMANCE COMPARISON (models vs ensemble)
    % =====================================================================

    if exist('ensemble_results.mat', 'file') && exist('forecast_metrics_summary.mat', 'file')

        fprintf('Plotting ensemble performance comparison...\n\n');

        E = load('ensemble_results.mat', 'results_ensemble');
        results_ensemble = E.results_ensemble;

        S = load('forecast_metrics_summary.mat', 'results');
        results = S.results;

        fig = figure('Name', 'Ensemble Performance', ...
                     'Position', [100, 100, 1000, 400], 'Visible', 'off');

        metrics = {'MAE', 'MSE', 'Coverage', 'WIS'};

        bar_labels = [models, {'Ensemble'}];

        for i = 1:numel(metrics)
            subplot(1, numel(metrics), i);

            vals = nan(K + 1, 1);

            for m = 1:K
                idx = strcmp(results.Model, models{m});
                vals(m) = mean(results.(sprintf('%s_forecast', metrics{i}))(idx), 'omitnan');
            end

            vals(K + 1) = mean(results_ensemble.(metrics{i}), 'omitnan');

            b = bar(vals);
            b.FaceColor = 'flat';
            for m = 1:K
                b.CData(m, :) = model_colors{m};
            end
            b.CData(K + 1, :) = ens_color;

            set(gca, 'XTick', 1:(K + 1), 'XTickLabel', bar_labels);
            xtickangle(30);
            ylabel(metrics{i});
            title(metrics{i});

            set(gca, 'FontSize', 10);
            grid on;
        end

        sgtitle('Performance Comparison: Individual Models vs Linear Pool Ensemble');
        save_figure(fig, fullfile(config.plots_folder, 'ensemble_comparison'));

    else
        fprintf('Ensemble or metrics results not found.\n');
        fprintf('Run collect_metrics_simple() and create_ensemble_simple() first.\n\n');
    end

    % =====================================================================
    % 5. SUMMARY METRICS BY WINDOW
    % =====================================================================

    fprintf('Plotting metrics by window...\n\n');

    if exist('forecast_metrics_summary.mat', 'file')

        S = load('forecast_metrics_summary.mat', 'results');
        results = S.results;

        have_ens = exist('ensemble_results.mat', 'file');
        if have_ens
            E = load('ensemble_results.mat', 'results_ensemble');
            results_ensemble = E.results_ensemble;
        end

        fig = figure('Name', 'Forecast Metrics by Window', ...
                     'Position', [100, 100, 1000, 800], 'Visible', 'off');

        metrics = {'MAE_forecast', 'MSE_forecast', 'Coverage_forecast', 'WIS_forecast'};
        ens_metrics = {'MAE', 'MSE', 'Coverage', 'WIS'};
        metric_labels = {'MAE', 'MSE', 'Coverage (%)', 'WIS'};

        all_markers = {'o-', 's-', 'd-', '^-', 'v-', 'p-', 'h-', '*-'};

        for i = 1:numel(metrics)
            subplot(2, 2, i);
            hold on;

            for m = 1:K
                idx = strcmp(results.Model, models{m});
                mk = all_markers{mod(m - 1, numel(all_markers)) + 1};
                plot(results.Window(idx), results.(metrics{i})(idx), mk, ...
                     'Color', model_colors{m}, 'LineWidth', 2, ...
                     'MarkerSize', 8, 'MarkerFaceColor', model_colors{m});
            end

            if have_ens
                plot(results_ensemble.Window, results_ensemble.(ens_metrics{i}), '-', ...
                     'Color', ens_color, 'LineWidth', 3, 'Marker', 'o', ...
                     'MarkerSize', 9, 'MarkerFaceColor', ens_color);
            end

            xlabel('Window');
            ylabel(metric_labels{i});
            title(metric_labels{i});

            if have_ens
                legend([models, {'Ensemble'}], 'Location', 'best');
            else
                legend(models, 'Location', 'best');
            end

            set(gca, 'FontSize', 10);
            grid on;
            box on;
        end

        sgtitle('Forecast Performance by Window');
        save_figure(fig, fullfile(config.plots_folder, 'metrics_by_window'));

    end

    fprintf('\n========================================================\n');
    fprintf('   PLOTTING COMPLETE!\n');
    fprintf('========================================================\n');
    fprintf('   All plots saved to: %s/\n', config.plots_folder);
    fprintf('========================================================\n\n');

end

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

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
    grid on;
    box on;
end

function path = ensemble_csv_path(tstart, config)
    path = fullfile(config.output_folder, sprintf( ...
        'Forecast-Ensemble-%s-tstart-%d-calibrationperiod-%d-horizon-%d.csv', ...
        config.ens_label, tstart, config.calib_period, config.forecast_horizon));
end

function csv_file = find_forecast_csv(model_name, tstart, config)
    % Find the forecast CSV for one model and window.
    %
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

function save_figure(fig, filepath)
    % Save figure as both PNG and FIG, then close

    saveas(fig, [filepath '.png']);
    saveas(fig, [filepath '.fig']);
    close(fig);

    fprintf('    Saved: %s.png and %s.fig\n', filepath, filepath);
end

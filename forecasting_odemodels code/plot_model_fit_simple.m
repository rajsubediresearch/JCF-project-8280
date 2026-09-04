% =========================================================================
% PLOT MODEL FIT WITH PARAMETER HISTOGRAMS
% For Teaching/Demonstration Purposes
% =========================================================================
% This script creates plots showing:
%   - Top row: Parameter histograms with 95% CI
%   - Bottom: Model fit (calibration) + forecast with different colored PIs
%
% =========================================================================

function plot_model_fit_simple()
    
    fprintf('\n========================================================\n');
    fprintf('   PLOTTING MODEL FIT WITH PARAMETERS\n');
    fprintf('========================================================\n\n');
    
    % =====================================================================
    % CONFIGURATION
    % =====================================================================
    
    config.dataset_file = 'ebola-scenario-1';
    config.dataset_name = 'Ebola Scenario 1';
    config.calib_period = 20;
    config.forecast_horizon = 4;
    config.n_windows = 5;
    config.window_shift = 4;
    config.output_folder = './output';
    config.plots_folder = './plots';
    
    % Create plots folder if it doesn't exist
    if ~exist(config.plots_folder, 'dir')
        mkdir(config.plots_folder);
        fprintf('Created folder: %s\n\n', config.plots_folder);
    end
    
    % Models (parameter names will be read from CSV)
    models = {'GLM', 'Richards', 'Gompertz'};
    
    % Colors
    calib_color = [1, 0.6, 0.6];     % Light red/salmon for calibration
    forecast_color = [0.6, 0.8, 1];  % Light blue for forecast
    
    % =====================================================================
    % LOAD DATA
    % =====================================================================
    
    fprintf('Loading data...\n');
    data = load(['./input/' config.dataset_file '.txt']);
    fprintf('  Loaded %d time points\n\n', size(data, 1));
    
    % =====================================================================
    % PROCESS EACH MODEL AND WINDOW
    % =====================================================================
    
    for m = 1:length(models)
        
        model_name = models{m};
        model_name_file = [model_name ' model'];
        
        fprintf('Processing %s model...\n', model_name);
        
        % ---------------------------------------------------------
        % Read parameter names from rolling window CSV
        % ---------------------------------------------------------
        
        param_csv_pattern = sprintf('parameters-rollingwindow-model_name-%s-*.csv', model_name_file);
        param_csv_files = dir(fullfile(config.output_folder, param_csv_pattern));
        
        if isempty(param_csv_files)
            fprintf('  Parameter CSV not found, skipping model...\n\n');
            continue;
        end
        
        param_csv = readtable(fullfile(config.output_folder, param_csv_files(1).name), ...
                              'VariableNamingRule', 'preserve');
        col_names = param_csv.Properties.VariableNames;
        
        % Extract parameter names (exclude 'time', CI bounds, and X0)
        params = {};
        for c = 2:length(col_names)
            name = col_names{c};
            % Skip CI bounds and X0
            if ~contains(name, '95%CI') && ~contains(name, 'X0')
                params{end+1} = name;
            end
        end
        
        n_params = length(params);
        fprintf('  Parameters: %s\n', strjoin(params, ', '));
        
        for w = 1:config.n_windows
            
            tstart = 1 + (w - 1) * config.window_shift;
            
            fprintf('  Window %d (tstart = %d)\n', w, tstart);
            
            % ---------------------------------------------------------
            % Load bootstrap parameter samples from MAT file
            % ---------------------------------------------------------
            
            mat_pattern = sprintf('Forecast-ODEModel-%s.txt-model_name-%s*tstart-%d*calibrationperiod-%d*.mat', ...
                config.dataset_file, model_name_file, tstart, config.calib_period);
            mat_files = dir(fullfile(config.output_folder, mat_pattern));
            
            if isempty(mat_files)
                fprintf('    MAT file not found, skipping...\n');
                continue;
            end
            
            mat_data = load(fullfile(config.output_folder, mat_files(1).name));
            
            % Get bootstrap parameter samples (Phatss_model1)
            if ~isfield(mat_data, 'Phatss_model1')
                fprintf('    No bootstrap samples found, skipping...\n');
                continue;
            end
            
            Phatss = mat_data.Phatss_model1;  % M x n_params matrix
            
            % Get forecast trajectories
            if ~isfield(mat_data, 'forecast_model12')
                fprintf('    No forecast trajectories found, skipping...\n');
                continue;
            end
            
            forecast_curves = mat_data.forecast_model12;  % T x M matrix
            
            % ---------------------------------------------------------
            % Create figure
            % ---------------------------------------------------------
            
            fig = figure('Name', sprintf('%s Window %d', model_name, w), ...
                         'Position', [100, 100, 900, 700], 'Visible', 'off');
            
            % ---------------------------------------------------------
            % Top row: Parameter histograms (skip fixed parameters)
            % ---------------------------------------------------------
            
            % First, identify which parameters have variation
            params_to_plot = {};
            param_indices = [];
            for p = 1:n_params
                param_samples = Phatss(:, p);
                param_LB = quantile(param_samples, 0.025);
                param_UB = quantile(param_samples, 0.975);
                
                % Only plot if parameter has variation (not fixed)
                if abs(param_UB - param_LB) > 1e-6
                    params_to_plot{end+1} = params{p};
                    param_indices(end+1) = p;
                end
            end
            
            n_params_plot = length(params_to_plot);
            
            if n_params_plot == 0
                fprintf('    No variable parameters to plot histograms\n');
                n_params_plot = 1;  % Avoid division by zero
            end
            
            for pp = 1:length(params_to_plot)
                subplot(2, n_params_plot, pp);
                
                p = param_indices(pp);
                param_samples = Phatss(:, p);
                
                % Calculate statistics
                param_median = median(param_samples);
                param_LB = quantile(param_samples, 0.025);
                param_UB = quantile(param_samples, 0.975);
                
                % Plot histogram
                histogram(param_samples, 20, 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
                
                xlabel(params_to_plot{pp});
                ylabel('Frequency');
                
                % Title with CI
                if param_UB > 1000
                    title(sprintf('%s (CI: %.2e, %.2e)', params_to_plot{pp}, param_LB, param_UB));
                else
                    title(sprintf('%s (CI: %.3f, %.3f)', params_to_plot{pp}, param_LB, param_UB));
                end
                
                set(gca, 'FontSize', 9);
                grid on;
                box on;
            end
            
            % ---------------------------------------------------------
            % Bottom: Model fit + forecast (spans all columns)
            % ---------------------------------------------------------
            
            subplot(2, max(n_params_plot, 1), (max(n_params_plot, 1) + 1):(2 * max(n_params_plot, 1)));
            hold on;
            
            % Time vectors (relative to window start)
            t_total = config.calib_period + config.forecast_horizon;
            t_rel = 1:t_total;  % Relative time
            t_calib = 1:config.calib_period;
            t_forecast = (config.calib_period + 1):t_total;
            
            % Get observed data
            obs_start = tstart;
            obs_end = min(tstart + t_total - 1, size(data, 1));
            obs_data = data(obs_start:obs_end, 2);
            t_obs = 1:length(obs_data);
            
            % Calculate median and 95% PI from bootstrap curves
            median_curve = median(forecast_curves(1:t_total, :), 2);
            LB_curve = quantile(forecast_curves(1:t_total, :)', 0.025)';
            UB_curve = quantile(forecast_curves(1:t_total, :)', 0.975)';
            
            % Ensure non-negative
            LB_curve = max(LB_curve, 0);
            UB_curve = max(UB_curve, 0);
            
            % Plot calibration period PI (salmon/red)
            fill([t_calib, fliplr(t_calib)], ...
                 [LB_curve(t_calib)', fliplr(UB_curve(t_calib)')], ...
                 calib_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            
            % Plot forecast period PI (blue)
            fill([t_forecast, fliplr(t_forecast)], ...
                 [LB_curve(t_forecast)', fliplr(UB_curve(t_forecast)')], ...
                 forecast_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            
            % Plot median curve
            plot(t_rel, median_curve, 'k-', 'LineWidth', 2);
            
            % Plot observed data
            % Calibration points (circles)
            calib_idx = t_obs <= config.calib_period;
            plot(t_obs(calib_idx), obs_data(calib_idx), 'ko', ...
                 'MarkerSize', 8, 'MarkerFaceColor', 'k');
            
            % Forecast points (triangles)
            forecast_idx = t_obs > config.calib_period;
            if any(forecast_idx)
                plot(t_obs(forecast_idx), obs_data(forecast_idx), 'k^', ...
                     'MarkerSize', 8, 'MarkerFaceColor', 'k');
            end
            
            % Vertical line at calibration end
            xline(config.calib_period, 'k--', 'LineWidth', 1.5);
            
            % Labels
            xlabel('Time (relative to window start)');
            ylabel('Cases');
            title(sprintf('Window %d - %s fit (Calibration + Forecast)', w, model_name));
            
            % Set axis limits
            xlim([0, t_total + 1]);
            ylim([0, max(UB_curve) * 1.1]);
            
            set(gca, 'FontSize', 11);
            grid on;
            box on;
            
            % ---------------------------------------------------------
            % Save figure
            % ---------------------------------------------------------
            
            save_figure(fig, fullfile(config.plots_folder, ...
                sprintf('%s_window%d_fit', lower(model_name), w)));
            
        end
        
        fprintf('\n');
    end
    
    fprintf('========================================================\n');
    fprintf('   PLOTTING COMPLETE!\n');
    fprintf('========================================================\n');
    fprintf('   All plots saved to: %s/\n', config.plots_folder);
    fprintf('========================================================\n\n');
    
end

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function save_figure(fig, filepath)
    % Save figure as both PNG and FIG, then close
    
    saveas(fig, [filepath '.png']);
    saveas(fig, [filepath '.fig']);
    close(fig);
    
    fprintf('    Saved: %s.png and %s.fig\n', filepath, filepath);
end
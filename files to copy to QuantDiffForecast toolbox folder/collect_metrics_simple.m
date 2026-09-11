% =========================================================================
% COLLECT FORECAST METRICS - SIMPLIFIED VERSION
% For Teaching/Demonstration Purposes
% =========================================================================
% This script collects calibration and forecast performance metrics
% (MAE, MSE, Coverage, WIS) and model selection criteria (AICc) from
% pre-computed CSV output files.
%
% Author: Raj Subedi
% =========================================================================

function collect_metrics_simple()
    
    fprintf('\n========================================================\n');
    fprintf('    METRICS COLLECTION\n');
    fprintf('========================================================\n\n');
    
    % =====================================================================
    % CONFIGURATION - Modify these for your dataset
    % =====================================================================
    
    config.dataset_name = 'Ebola Scenario 1';  % As appears in filenames
    config.calib_period = 20;                  % Calibration window (weeks)
    config.forecast_horizon = 4;               % Forecast horizon (weeks)
    config.n_windows = 5;                      % Number of forecast windows
    config.window_shift = 4;                   % Weeks between window starts
    config.output_folder = './output';         % Where CSV files are stored
    
    % Models to evaluate
    models = {'GLM', 'Richards', 'Gompertz'};
    
    fprintf('Configuration:\n');
    fprintf('  Dataset: %s\n', config.dataset_name);
    fprintf('  Calibration period: %d weeks\n', config.calib_period);
    fprintf('  Forecast horizon: %d weeks\n', config.forecast_horizon);
    fprintf('  Windows: %d (shift = %d weeks)\n\n', config.n_windows, config.window_shift);
    
    % =====================================================================
    % COLLECT METRICS FOR EACH MODEL AND WINDOW
    % =====================================================================
    
    % Initialize results table
    results = table();
    
    for w = 1:config.n_windows
        
        % Calculate window start time
        tstart = 1 + (w - 1) * config.window_shift;
        
        fprintf('Window %d (tstart = %d):\n', w, tstart);
        
        for m = 1:length(models)
            model = models{m};
            
            % Initialize row for this model-window combination
            row = struct();
            row.Window = w;
            row.tstart = tstart;
            row.Model = {model};  % Cell array for table compatibility
            
            % ---------------------------------------------------------
            % 1. Read CALIBRATION metrics
            % ---------------------------------------------------------
            calib_file = find_csv('calibration', model, tstart, config);
            
            if ~isempty(calib_file)
                calib = readtable(calib_file, 'VariableNamingRule', 'preserve');
                row.MAE_calib = calib.MAE(1);
                row.MSE_calib = calib.MSE(1);
                row.Coverage_calib = calib.("Coverage 95%PI")(1);
                row.WIS_calib = calib.WIS(1);
            else
                warning('  Calibration file not found for %s', model);
                row.MAE_calib = NaN;
                row.MSE_calib = NaN;
                row.Coverage_calib = NaN;
                row.WIS_calib = NaN;
            end
            
            % ---------------------------------------------------------
            % 2. Read FORECAST metrics
            % ---------------------------------------------------------
            forecast_file = find_csv('forecasting', model, tstart, config);
            
            if ~isempty(forecast_file)
                forecast = readtable(forecast_file, 'VariableNamingRule', 'preserve');
                % Average across all forecast horizons
                row.MAE_forecast = mean(forecast.MAE);
                row.MSE_forecast = mean(forecast.MSE);
                row.Coverage_forecast = mean(forecast.("Coverage 95%PI"));
                row.WIS_forecast = mean(forecast.WIS);
            else
                warning('  Forecast file not found for %s', model);
                row.MAE_forecast = NaN;
                row.MSE_forecast = NaN;
                row.Coverage_forecast = NaN;
                row.WIS_forecast = NaN;
            end
            
            % ---------------------------------------------------------
            % 3. Read AICc (model selection criterion)
            % ---------------------------------------------------------
            aic_file = find_aic_csv(model, tstart, config);
            
            if ~isempty(aic_file)
                aic = readtable(aic_file, 'VariableNamingRule', 'preserve');
                row.AICc = aic.AICc(1);
            else
                warning('  AICc file not found for %s', model);
                row.AICc = NaN;
            end
            
            % Add row to results
            results = [results; struct2table(row)];
            
            fprintf('  %s: MAE=%.2f, MSE=%.2f, Coverage=%.1f%%, WIS=%.2f, AICc=%.2f\n', ...
                model, row.MAE_forecast, row.MSE_forecast, row.Coverage_forecast, ...
                row.WIS_forecast, row.AICc);
        end
        fprintf('\n');
    end
    
    % =====================================================================
    % SUMMARY STATISTICS
    % =====================================================================
    
    fprintf('========================================================\n');
    fprintf('   SUMMARY BY MODEL (averaged across windows)\n');
    fprintf('========================================================\n\n');
    
    fprintf('%-10s | MAE_fore | MSE_fore | Coverage | WIS_fore | AICc\n', 'Model');
    fprintf('-----------|----------|----------|----------|----------|----------\n');
    
    for m = 1:length(models)
        idx = strcmp(results.Model, models{m});
        fprintf('%-10s | %8.2f | %8.2f | %7.1f%% | %8.2f | %8.2f\n', ...
            models{m}, ...
            mean(results.MAE_forecast(idx), 'omitnan'), ...
            mean(results.MSE_forecast(idx), 'omitnan'), ...
            mean(results.Coverage_forecast(idx), 'omitnan'), ...
            mean(results.WIS_forecast(idx), 'omitnan'), ...
            mean(results.AICc(idx), 'omitnan'));
    end
    
    fprintf('\n');
    
    % =====================================================================
    % SAVE OUTPUT
    % =====================================================================
    
    % Save as CSV
    output_csv = 'forecast_metrics_summary.csv';
    writetable(results, output_csv);
    fprintf('Saved: %s\n', output_csv);
    
    % Save as MAT file
    output_mat = 'forecast_metrics_summary.mat';
    save(output_mat, 'results', 'config', 'models');
    fprintf('Saved: %s\n\n', output_mat);
    
    % Display table preview
    fprintf('========================================================\n');
    fprintf('   RESULTS TABLE PREVIEW\n');
    fprintf('========================================================\n\n');
    disp(results);
    
    fprintf('Done!\n\n');
end

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function filepath = find_csv(type, model, tstart, config)
    % Find calibration or forecast performance CSV file
    %
    % type: 'calibration' or 'forecasting'
    % model: 'GLM', 'Richards', or 'Gompertz'
    
    pattern = sprintf('performance-%s-*model_name-%s*tstart-%d*calibrationperiod-%d*.csv', ...
        type, model, tstart, config.calib_period);
    
    files = dir(fullfile(config.output_folder, pattern));
    
    if isempty(files)
        filepath = '';
    else
        filepath = fullfile(config.output_folder, files(1).name);
    end
end

function filepath = find_aic_csv(model, tstart, config)
    % Find AICc CSV file
    
    pattern = sprintf('AICc-model_name-%s*tstart-%d*calibrationperiod-%d*.csv', ...
        model, tstart, config.calib_period);
    
    files = dir(fullfile(config.output_folder, pattern));
    
    if isempty(files)
        filepath = '';
    else
        filepath = fullfile(config.output_folder, files(1).name);
    end
end
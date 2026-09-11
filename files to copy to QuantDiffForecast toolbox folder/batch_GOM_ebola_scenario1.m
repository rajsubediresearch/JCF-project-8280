% EBOLA SCENARIO 1 - FIXED WINDOWS (GOMPERTZ MODEL)
% Fixed calibration period: 20 weeks %changed for this project
% Shift by forecast horizon: 4 weeks (non-overlapping forecasts)

fprintf('========================================\n');
fprintf('   EBOLA SCENARIO 1 - FIXED WINDOWS\n');
fprintf('   Model: GOMPERTZ\n');
fprintf('========================================\n\n');

% Configuration
calib_size = 20; %changed for this project
forecast_horizon = 4;
shift = 4;  % Shift by forecast horizon

% Load data
data = load('./input/ebola-scenario-1.txt');
max_time = size(data, 1);
n_windows = floor((max_time - calib_size - forecast_horizon + 1) / shift);

fprintf('Configuration:\n');
fprintf('  Calibration period: %d weeks (FIXED)\n', calib_size);
fprintf('  Forecast horizon: %d weeks\n', forecast_horizon);
fprintf('  Shift: %d weeks\n', shift);
fprintf('  Total data points: %d weeks\n', max_time);
fprintf('  Number of windows: %d\n\n', n_windows);

% Run forecasts
for window = 1:n_windows
    tstart = 1 + (window - 1) * shift;  % Start: 1, 5, 9, 13, 17, 21
    
    fprintf('[Window %2d/%2d] tstart=%2d, calib=weeks %2d-%2d, forecast=weeks %2d-%2d\n', ...
            window, n_windows, tstart, ...
            tstart, tstart + calib_size - 1, ...
            tstart + calib_size, tstart + calib_size + forecast_horizon - 1);
    
    Run_Forecasting_ODEModel(@options_forecast_GOM_ebola_scenario1, tstart, tstart, calib_size, forecast_horizon);
end

fprintf('\n========================================\n');
fprintf('✓ All GOMPERTZ fixed window runs completed\n');
fprintf('✓ Output files saved to ./output/\n');
fprintf('========================================\n');

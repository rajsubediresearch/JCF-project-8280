function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B, model, params,vars,getperformance, forecastingperiod,windowsize1,tstart1,tend1,printscreen1]=options_forecast
% Options file for Richards model
% Weekly data, fixed rolling windows

global method1

% Dataset
cadfilename1='JALISCO_2025-08-18_2026-06-29-trimmed';
caddisease='measles';
datatype='cases';

% Parameter estimation
method1=3;
dist1=3;

switch method1
    case 1
        dist1=1;
    case 3
        dist1=3;
    case 4
        dist1=4;
    case 5
        dist1=5;
end

numstartpoints=100;
B=300;

% Model: Richards
model.fc = @RICH;  % Function name in toolbox
model.name = 'Richards model';  % Display name for files/plots
params.label = {'r','a','K'};
params.LB = [0.05, 0.1, 1000];     % Lower bounds
params.UB = [5, 1, 10000];            % Upper bounds (adjusted for Ebola)
params.initial = [0.5, 0.5, 7500];    % Initial guesses (adjusted for Ebola)
params.fixed = [0, 0, 0];
params.fixI0 = 1;
params.composite='';
params.extra0='';

vars.label={'C'};
vars.initial=1;  %first case in series
vars.fit_index=1;
vars.fit_diff=1;

% Forecasting
getperformance=1;
forecastingperiod=3;

% Rolling window, will be overridden during execution
windowsize1=40;
tstart1=1;
tend1=1;
printscreen1=1;

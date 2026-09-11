function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B, model, params,vars,getperformance, forecastingperiod,windowsize1,tstart1,tend1,printscreen1]=options_forecast
% Options file for Richards model - Ebola Scenario 1
% Weekly data, fixed rolling windows

global method1

% Dataset
cadfilename1='ebola-scenario-1';
caddisease='Ebola Scenario 1';
datatype='cases';

% Parameter estimation
method1=0;
dist1=0;

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

numstartpoints=20;
B=300;

% Model: Richards
model.fc = @RICH;  % Function name in toolbox
model.name = 'Richards model';  % Display name for files
params.label = {'r','a','K'};
params.LB = [0.001, 0.05, 200];     % Lower bounds
params.UB = [1, 3, 50000];            % Upper bounds (adjusted for Ebola)
params.initial = [0.1, 1, 10000];    % Initial guesses (adjusted for Ebola)
params.fixed = [0, 0, 0];
params.fixI0 = 1;
params.composite='';
params.extra0='';

vars.label={'C'};
vars.initial=8;  % Ebola Scenario 1: first case count is 8
vars.fit_index=1;
vars.fit_diff=1;

% Forecasting
getperformance=1;
forecastingperiod=4;

% Rolling window
windowsize1=20; %edited for 20 weeks cal period for each window. windows will be altered by the batch file
tstart1=1;
tend1=1;
printscreen1=0;

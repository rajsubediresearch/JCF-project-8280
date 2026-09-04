function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B, model, params,vars,getperformance, forecastingperiod,windowsize1,tstart1,tend1,printscreen1]=options_forecast
% Options file
% Weekly data, fixed rolling windows

global method1

% Dataset
cadfilename1='JALISCO_2025-08-18_2026-06-29-trimmed';
caddisease='measles';
datatype='cases';

% Parameter estimation method
method1=3;  %NegBin
dist1=3;    %NegBin

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

% Model: GLM
model.fc = @GLM;
model.name = 'GLM model';
params.label = {'r','p','K'};
params.LB = [0.05, 0.1, 1000];      % Lower bounds
params.UB = [5, 1, 10000];           % Upper bounds
params.initial = [0.5, 0.5, 7500];   % Initial guesses
params.fixed = [0, 0, 0];
params.fixI0 = 1;
params.composite='';
params.extra0='';

vars.label={'C'};
vars.initial=1;  % first case count in the series
vars.fit_index=1;
vars.fit_diff=1;

% Forecasting
getperformance=1;
forecastingperiod=3;  % 3 weeks ahead

% Rolling window, will be overridden during execution with supplied values
windowsize1=40;
tstart1=1;
tend1=1;
printscreen1=1;

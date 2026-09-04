function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B, model, params,vars,getperformance, forecastingperiod,windowsize1,tstart1,tend1,printscreen1]=options_forecast
% Options file for GLM model - Ebola Scenario 1
% Weekly data, fixed rolling windows

global method1

% Dataset
cadfilename1='ebola-scenario-1';
caddisease='Ebola Scenario 1';
datatype='cases';

% Parameter estimation method
method1=0;  % 0=LSQ, 3=MLE NegBin
dist1=0;    % 0=Normal, 3=NegBin

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

% Model: GLM
model.fc = @GLM;
model.name = 'GLM model';
params.label = {'r','p','K'};
params.LB = [0.001, 0.1, 100];      % Lower bounds
params.UB = [2, 1, 15000];           % Upper bounds (adjusted for Ebola)
params.initial = [0.1, 0.5, 2000];   % Initial guesses (adjusted for Ebola)
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
forecastingperiod=4;  % 4 weeks ahead

% Rolling window (will be overridden by batch script)
windowsize1=20; %edited based on cal length of 20 weeks
tstart1=1;
tend1=1;
printscreen1=0;

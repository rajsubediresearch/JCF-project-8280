% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function dx=LM(t,x,params0,extra0)

% Logistic growth model (GrowthPredict flag1=3)
% parameters in order: r, K

dx=zeros(1,1);

dx(1,1)=params0(1)*x(1,1).*(1-(x(1,1)/params0(2)));

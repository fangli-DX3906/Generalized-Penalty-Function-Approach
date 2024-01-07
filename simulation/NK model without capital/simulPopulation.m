%**************************************************************
% MAIN_PROG - Simulating the Population
%
%**************************************************************
cd '/Users/fangli/Dropbox/Marco_Fang_project/SimpleNKGPFA/NK model without capital'
addpath '/Users/fangli/Dropbox/Marco_Fang_project/SimpleNKGPFA/NK model without capital/funcs'

clc
clear

%% 
%%
% run solveNK.m
run solveNK_HF.m
% run solveNK_Calvo.m
% run solveNK_Nonsep.m

close all

%% parameter
%%
ecolength = 100000;

%% simulation
%%
trueShocks = zeros(ecolength+1, 2);
ds = randn(1,1);
ss = randn(1,1);
trueShocks(1,:) = [ds, ss];
x0 = eta * [ss; 0; 0; ds];     % t=1 states = t=1 shocks
pd = length(x0);
MX = [gx; eye(pd)];
x= x0;
for t=1:ecolength
    Y(t, :) = (MX*x)';        % t=1 controls, the start of the series passed to gpfa
    ds = randn(1,1);
    ss = randn(1,1);
    trueShocks(t+1, :) = [ds, ss];  % t=2 random shocks
    x = hx*x+eta*[ss; 0; 0; ds];     % t=2 states
end
trueShocks = trueShocks(1:end-1, :);

% save './data/population.mat' Y
% save './data/trueShockSeries.mat' trueShocks
save './data/HF/population.mat' Y
save './data/HF/trueShockSeries.mat' trueShocks
% save './data/Calvo/population.mat' Y
% save './data/Calvo/trueShockSeries.mat' trueShocks
% save './data/Nonsep/population.mat' Y
% save './data/Nonsep/trueShockSeries.mat' trueShocks
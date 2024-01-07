%**************************************************************
% MAIN_PROG - Simulating 2000 short economies
%
%**************************************************************
%cd '/Users/fangli/Dropbox/Marco_Fang_project/SimpleNKGPFA/NK model without capital'
cd 'C:\Users\brianti\Dropbox\Share_Marco_Fang\SimpleNKGPFA\NK model without capital'
%addpath '/Users/fangli/Dropbox/Marco_Fang_project/SimpleNKGPFA/NK model without capital/funcs'
addpath 'C:\Users\brianti\Dropbox\Share_Marco_Fang\SimpleNKGPFA\NK model without capital\funcs'

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
nsim = 10;
ecolength = 24000;

%% simulation
%%
Ysim = cell(nsim,1);
Shocksim = cell(nsim, 1);

for isim = 1:nsim
    shockr = zeros(ecolength+1, 2);
    ds = randn(1,1);
    ss = randn(1,1);
    shockr(1, :) = [ds, ss];
    x0 = eta * [ss; 0; 0; ds];
    pd = length(x0);
    MX = [gx; eye(pd)];
    x= x0;
    
    for t=1:ecolength
        Y(t, :) = (MX*x)';
        ds = randn(1,1);
        ss = randn(1,1);
        shockr(t+1, :) = [ds, ss];
        x = hx*x+eta*[ss; 0; 0; ds];
    end
    
    Ysim{isim} = Y;
    Shocksim{isim} = shockr(1:end-1, :);
end

% save './data/economies2000.mat' Ysim
% save './data/shocksim2000.mat' Shocksim
save './data/HF/economies2000.mat' Ysim
save './data/HF/shocksim2000.mat' Shocksim
% save './data/Nonsep/economies2000.mat' Ysim
% save './data/Nonsep/shocksim2000.mat' Shocksim
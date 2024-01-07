%% Demand vs Supply shock
%%

clear
close all
clc
cd '/Users/fangli/Library/Mobile Documents/com~apple~CloudDocs/Empirics/thisOne'

%% params value                           
system_names = {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP', 'Employee'}; 
START = 1960;  % 1960-2021.9.30
END = 2020;
iscon = 1; 
istr = 0; 
IRhoriz = 51; 
lags = 1:8;

%% import the data
[dataset, txt, ~] = xlsread('2021data.xlsx');
tf = isreal(dataset);
if tf == 0
      fprintf('\n')
      warning('Dataset has complex variables in it.')
      dataset = real(dataset);
      fprintf('\n')
end
time_start  = dataset(1,1);
time_end = dataset(end,1);
Time = dataset(:,1);

% Assess names to each variable as an array
for i = 1:size(dataset,2)
      eval([txt{i} ' = dataset(:,i);']);
end

% Standard Transformation
Consumption = log(RealConsumptionNonDurables + RealConsumptionService);
Investment = log(RealConsumptionDurables + RealInvestment);
GDP = log(RealGDP);
CPI = log(CPI);
% Def = log(GDPDeflator);
TotalHours = log(TotalHours);
Employee = log(Employee); 
TFP = cumsum(TFP,'omitnan')/400;

% Build System (data set)
for i = 1:length(system_names)
      Y(:,i) = eval(system_names{i});
end

% Truncate data in case of NaN
[XXX, ~, ~] = truncate_data([Y Time]);
Y = XXX(:,1:end-1);
TimeXXX = XXX(:,end);
locSTART = find(TimeXXX == START);
locEND = find(TimeXXX == END);
TimeXXX = TimeXXX(1:locEND);
Y = Y(1:locEND, :);

if isempty(locSTART) == 1
      Y = Y(1:end,:);
      TimeXXX = TimeXXX(1:end);
      fprintf('\n')
      warning(['First available observation is in ',num2str(TimeXXX(1))])
      fprintf('\n')
else
      Y = Y(locSTART:end,:);
      TimeXXX = TimeXXX(locSTART:end);
end

%% 0. basic VAR estimation
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test the validity
Y0_ = Y(:, 1:end-1);
Y0_ = Y(:, [1 3]);
Y0_ = diff(Y0_);
% Y0_ = diff(Y0_);
system_names0_ = {'GDP', 'Investment'}; 
for i=1:size(lags,2)
    [~, ~, ~, AA0, ~, ~, ~] = estimVARCoeff(Y0_, lags(i), iscon, istr, system_names0_); 
    flag = test_stationarity(AA0'); 
    if flag == 0
        disp([num2str(i), 'lag case is stationary'])
    else
        disp([num2str(i), 'lag case is not stationary'])
    end
    clear temp_mat
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y0 = Y(:, 1:end-1);
system_names0 = {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'};
for i=1:size(lags,2)
    [~, ~, ~, AA0, ~, ~, ~] = estimVARCoeff(Y0, lags(i), iscon, istr, system_names0); 
    flag = test_stationarity(AA0');  % a bug in the sign, pay attention to the constant
    if flag == 0
        disp([num2str(i), 'lag case is stationary'])
    else
        disp([num2str(i), 'lag case is not stationary'])
    end
    clear temp_mat
end
% only lag = 1 is stationary

[Bcomp, ~, ~, Bpl_ev, VC_eps, Residmat, mVAR] = estimVARCoeff(Y0, 4, iscon, istr, system_names0); 
IRFs = calcIR(mVAR, Residmat, 'orthogonalized', IRhoriz, 1, "btsp", 0.95);
plotIR(mVAR, IRFs, {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}, [1,2,3,4,5,6,7], 1, {'1', '2', '3','4','5','6','7'});
fig = figure(1);
fig.Color = 'w'; 
fig.WindowState = 'maximized';
% A = mVAR.auxMats.cholVC;
% IRFs_tt = calcOrthoIRFs(mVAR, A, IRhoriz, NaN);
% IRFs_tt = mkIRFStruct(mVAR, IRFs_tt, NaN, NaN, 0.95);
% plotIR(mVAR, IRFs_tt, {'GDP', 'CPI'}, [1, 2], 1, {'shock 1', 'shock 2'});

% identification
% initDelta = 0.988;   % for DEN case
initDelta = 0;     % for no DEN case
[GAMMA, mVAR] = solver(mVAR, {{'GDP', 'CPI'}, {'GDP', 'CPI'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], initDelta, [1,2], 1, [], 1);  
A_mat = mVAR.auxMats.PFMats;
IRFs_PF = calcOrthoIRFs(mVAR, A_mat, IRhoriz, NaN);
plotIR(mVAR, IRFs_PF, {'GDP', 'CPI'}, [1, 2], 1, {'demand shock', 'supply shock'});

% bootstrap for CI
[uci, lci] = bootstrapCI(mVAR, 50, 0.95, IRhoriz, []);
% [uci, lci] = bootstrapCI_(mVAR, 100, 0.95, IRhoriz, []);
IRFs_ = mkIRFStruct(mVAR, IRFs_PF, uci, lci, 0.95);
plotIR(mVAR, IRFs_, {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}, [1, 2], 1, {'demand shock', 'supply shock'});

%% 1. Adding inflation to the system (replacing 'CPI')
%%
Inflation = log(Y(2:end, 2)) - log(Y(1:end-1, 2));
Y1 = Y(2:end, 1:end-1);
Y1(:, 2) = Inflation;
system_names1 = {'GDP', 'Inflation', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}; 
for i=1:size(lags,2)
    [~, ~, ~, AA1, ~, ~, ~] = estimVARCoeff(Y1, lags(i), iscon, istr, system_names1); 
    flag = test_stationarity(AA1');
    if flag == 0
        disp([num2str(i), 'lag case is stationary'])
    else
        disp([num2str(i), 'lag case is not stationary'])
    end
    clear temp_mat
end
% only lag = 1 is stationary

[Bcomp, ~, ~, Bpl_ev, VC_eps, Residmat, mVAR1] = estimVARCoeff(Y1, 1, iscon, istr, system_names1); 
IRFs1 = calcIR(mVAR1, Residmat, 'orthogonalized', IRhoriz, 1, "btsp", 0.95);
plotIR(mVAR1, IRFs1, {'GDP', 'Inflation'}, [1, 2], 1, {'shock 1', 'shock 2'});

initDelta = 0;     % for no DEN case
[GAMMA, mVAR1] = solver(mVAR1, {{'GDP', 'Inflation'}, {'GDP', 'Inflation'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], initDelta, [1,2], 1, [], 1);  
A_mat1 = mVAR1.auxMats.PFMats;
IRFs_PF1 = calcOrthoIRFs(mVAR1, A_mat1, IRhoriz, NaN);
plotIR(mVAR1, IRFs_PF1, {'GDP', 'Inflation'}, [1, 2], 1, {'demand shock', 'supply shock'});

[uci, lci] = bootstrapCI(mVAR1, 50, 0.95, IRhoriz, []);
IRFs_1 = mkIRFStruct(mVAR1, IRFs_PF1, uci, lci, 0.95);
plotIR(mVAR1, IRFs_1, {'GDP', 'Inflation', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}, [1, 2], 1, {'demand shock', 'supply shock'});

%% 2. linearly detrend all variables (not using inflation)
%%
timeIdx = 1:size(Y, 1);
XX = [ones(size(Y, 1), 1) timeIdx'];
Y2 = zeros(size(Y, 1), size(Y, 2)-1);
for i = 1: size(Y2, 2)
    [~, ~, Y2(:, i)] = regress(Y(:, i), XX);
end
system_names2 = {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}; 
for i=1:size(lags,2)
    [~, ~, ~, AA2, ~, ~, ~] = estimVARCoeff(Y2, lags(i), iscon, istr, system_names2); 
    flag = test_stationarity(AA2');
    if flag == 0
        disp([num2str(i), 'lag case is stationary'])
    else
        disp([num2str(i), 'lag case is not stationary'])
    end
    clear temp_mat
end
% no stationary case

[Bcomp, ~, ~, Bpl_ev, VC_eps, Residmat, mVAR2] = estimVARCoeff(Y2, 1, iscon, istr, system_names2); 
IRFs2 = calcIR(mVAR2, Residmat, 'orthogonalized', IRhoriz, 1, "btsp", 0.95);
plotIR(mVAR2, IRFs2, {'GDP', 'CPI'}, [1, 2], 1, {'shock 1', 'shock 2'});

initDelta = 0;     % for no DEN case
[GAMMA, mVAR2] = solver(mVAR2, {{'GDP', 'CPI'}, {'GDP', 'CPI'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], initDelta, [1,2], 1, [], 1);  
A_mat = mVAR2.auxMats.PFMats;
IRFs_PF2 = calcOrthoIRFs(mVAR2, A_mat, IRhoriz, NaN);
plotIR(mVAR2, IRFs_PF2, {'GDP', 'CPI'}, [1, 2], 1, {'demand shock', 'supply shock'});

[uci, lci] = bootstrapCI(mVAR2, 50, 0.95, IRhoriz, []);
IRFs_2 = mkIRFStruct(mVAR2, IRFs_PF2, uci, lci, 0.95);
plotIR(mVAR2,  IRFs_2, {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}, [1, 2], 1, {'demand shock', 'supply shock'});

%% 3.Removing population
%%
system_names3 = {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}; 
Y3 = Y;
Y3(:, 1) = Y3(:, 1) - Y3(:, end);
Y3(:, 5) = Y3(:, 5) - Y3(:, end);
Y3(:, 3) = Y3(:, 3) - Y3(:, end);
Y3(:, 4) = Y3(:, 4) - Y3(:, end);
Y3 = Y3(:, 1:end-1);
for i=1:size(lags,2)
    [~, ~, ~, AA3, ~, ~, ~] = estimVARCoeff(Y3, lags(i), iscon, istr, system_names3); 
    flag = test_stationarity(AA3');
    if flag == 0
        disp([num2str(i), 'lag case is stationary'])
    else
        disp([num2str(i), 'lag case is not stationary'])
    end
    clear temp_mat
end
% only lag = 1 is stationary

[Bcomp, ~, ~, Bpl_ev, VC_eps, Residmat, mVAR3] = estimVARCoeff(Y3, 1, iscon, istr, system_names3); 
IRFs3 = calcIR(mVAR3, Residmat, 'orthogonalized', IRhoriz, 1, "btsp", 0.95);
plotIR(mVAR3, IRFs3, {'GDP', 'CPI'}, [1, 2], 1, {'shock 1', 'shock 2'});

initDelta = 0;     % for no DEN case
[GAMMA, mVAR3] = solver(mVAR3, {{'GDP', 'CPI'}, {'GDP', 'CPI'}}, [[1,1]; [1,1]], [[1,1]; [1,-1]], initDelta, [1,2], 1, [], 1);  
A_mat = mVAR3.auxMats.PFMats;
IRFs_PF3 = calcOrthoIRFs(mVAR3, A_mat, IRhoriz, NaN);
plotIR(mVAR3, IRFs_PF3, {'GDP', 'CPI'}, [1, 2], 1, {'demand shock', 'supply shock'});

[uci, lci] = bootstrapCI(mVAR3, 50, 0.95, IRhoriz, []);
IRFs_3 = mkIRFStruct(mVAR3, IRFs_PF3, uci, lci, 0.95);
plotIR(mVAR3,  IRFs_3, {'GDP', 'CPI', 'Investment', 'Consumption', 'TotalHours', 'ShadowRate', 'TFP'}, [1, 2], 1, {'demand shock', 'supply shock'});
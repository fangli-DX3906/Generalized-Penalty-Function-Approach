%
%
% Read monthly data and aggregate them at a quarterly frequency
%
% Code by Brianti, M. 
% Today 06/24/2022
%

clear
close all
clc
path(pathdef);
addpath('C:\Users\brianti\Dropbox\Job_Market_Paper\Replication Brianti 2021\Empirics 2022\data');  % Windows


% Data
filename                    = 'quarterly_24june_2022'; %don-t have {'#N/A'} in the excel, just blank
sheet                       = 'monthly';
range                       = 'B2:O750';
range2                      = 'B1:O1'; % it can be done in one step but it was not working for Marco
[RawData,~]                 = xlsread(filename,sheet,range);
[~,var_names]               = xlsread(filename,sheet,range2); % it can be done in one step but it was not working for Marco

% Assess name to each variable
for iv = 1:size(RawData,2)
    eval([var_names{iv} ' = RawData(:,iv);']);
end

% Aggregate via simple average
iq = 0;
for im = 1:length(RawData)-2
    if round((im+2)/3) == (im+2)/3
        iq = iq + 1;
        RawDataq(iq,:) = mean(RawData(im:im+2,:),1);
        RawDataq(iq,1) = RawDataq(iq,1) - 1/12;
    end
end

% Add back to excel
filename = 'aggregated_data';
sheet    = 1;
xlswrite(filename,var_names,sheet,'A1');
xlswrite(filename,RawDataq,sheet,'A2');




function [flexp, setp] = paramset_HF() 

setp.bet              = 0.99;    % 0.99 Never call anything beta...it's a matlab function
setp.alp              = 0.33;    % 0.33 share of physical capital. Thus, 0.66 is the share of labor income in this model without physical capital
setp.sigmc          = 2;       % CRRA parameter of consumption between (1 2] 
setp.sigmn          = 2;       % CRRA parameter of leaisure between (1 2] 
setp.nss              = 0.33;    % share of hours worked in ss
setp. bb              = 0.9;
setp.siga             = 0.01;    % 0.01 variance of tech shock (keep one percent for now)
setp.sigt             = 0.01;    % 0.01 variance of government spending shock (keep one percent for now)
% setp.del             = 0.025;   % physical capital depreciation rate
setp.eta               = 2;       % elasticity of substitution across goods yi
setp.pp                = 5;     % Taylor coefficient from the Taylor rule (monetary policy rule)
setp.gamp           = 15;      % Price adjustment cost (a la Rotemberg) parameter
setp.rhoa            = 0.75;    % 0.95 persistence of tech shock
setp.rhot            = 0.9;    % persistence of government spending shock
setp.rhoc            =0.9;
setp.rhob            =0.75;
setp.rhoi            = 0.9;    % interest rate smoothing parameter
setp.psii             = 9999;    % Unused since psi is estimated in steady state in function of nss = 1/3

% Params to be estimated via IRF matching
flexp.name           = nan;     % empty for now

end



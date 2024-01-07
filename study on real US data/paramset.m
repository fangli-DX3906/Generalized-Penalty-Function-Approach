function [flexp, setp] = paramset(whichSet)

if whichSet ==1 

    setp.bet             = 0.99;                      % 0.99 Never call anything beta...it's a matlab function
    setp.alp             = 0.33;                      % 0.33 share of physical capital
    setp.sigmc         = 2;                          % CRRA parameter of consumption between (1 2] 
    setp.sigmn         = 2;                          % CRRA parameter of leaisure between (1 2] 
    setp.nss             = 0.33;                     % share of hours worked in ss
    setp.del             = 0.025;                    % physical capital depreciation rate
    setp.eta             = 2;                          % elasticity of substitution across goods yi
    setp.pp              = 1.5;                       % Taylor coefficient from the Taylor rule (monetary policy rule)
    setp.gamp         = 2;                        % Price adjustment cost (a la Rotemberg) parameter
    % param1, gamp=10, param0, gamp=1 
    setp.rhoa           = 0.95;                     % 0.95 persistence of tech shock
    setp.rhot            = 0.9;                      % persistence of government spending shock
    setp.rhoi            = 0.75;                     % interest rate smoothing parameter
    setp.rhoc           = 0.9;                       % persistence of the shock to monetary policy
    setp.rhob           = 0.75;                     % persistence of the shock to preference 
    setp.siga            = 0.01;                     % 0.01 variance of tech shock (keep one percent for now)
    setp.sigt             = 0.01;                     % 0.01 variance of government spending shock (keep one percent for now)
    setp.sigb            = 0.01;                     % 0.01 variance of preference shock (keep one percent for now)
    setp.sigc            = 0.01;                     % 0.01 variance of monetary policy shock (keep one percent for now)
    setp.psii            = 9999;                     % Unused since psi is estimated in steady state in function of nss = 1/3

else
    % new set of params
    setp.bet             = 0.99;    
    setp.alp             = 0.33;    
    setp.sigmc         = 2;      
    setp.sigmn         = 2;      
    setp.nss             = 0.33;   
    setp.siga            = 0.01;    
    setp.sigt            = 0.01;    
    setp.del             = 0.025;   
    setp.eta             = 2;       
    setp.pp              = 2;    
    setp.gamp         = 10;      
    setp.rhoa            = 0.75;   
    setp.rhot            = 0.90;   
    setp.rhoi            = 0;    
    setp.rhob            = 0.75;     
    setp.psii            = 9999;    
    setp.rhoc            = 0;       
    
end

% Params to be estimated via IRF matching
flexp.name           = nan;     % empty for now

end


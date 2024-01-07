function [flexp, setp] = paramset_Nonsep() 

setp.bet             = 0.99;   
setp.alp            = 0.33;
setp.sga             = 2;
setp.gga            = 0.5;
% setp.nss             = 0.33;    
setp.eta             = 2;       
setp.pp              = 2;     
setp.theta            = 0.75;      
setp.rhoa            = 0.75;    
setp.rhot            = 0.90;    
setp.rhoc            = 0.9;    
setp.rhob            = 0.75; 
setp.rhoi            = 0.9;    
setp.gamp           = 15;      % Price adjustment cost (a la Rotemberg) parameter
setp.psii              =9999;

% Params to be estimated via IRF matching
flexp.name           = nan;     % empty for now

end
function [flexp, setp] = paramset_Calvo() 

setp.bet             = 0.99;    
setp.sigmc           = 2;       
setp.sigmn           = 2;       
setp.nss             = 0.33;    
setp.eta             = 2;       
setp.pp              = 2;     
setp.theta            = 0.75;      
setp.rhoa            = 0.75;    
setp.rhot            = 0.90;    
setp.rhoc            = 0.9;    
setp.rhob            = 0.75; 
setp.rhoi            = 0.9;    
setp.psii            = 9999;    

% Params to be estimated via IRF matching
flexp.name           = nan;     % empty for now

end



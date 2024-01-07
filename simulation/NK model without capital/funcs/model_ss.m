function [ss, setp, flexp] = model_ss(flexp,setp)

% [ss, param] = model_ss(param)
% Return the steady state of the model (computed analytically)
% Code by Marco Brianti, University of Alberta
% January 1, 2022

% Parameters (flexp), from char to num with the same name. No flexp now.
nms = fieldnames(flexp);
for j = 1:length(nms)
    eval([nms{j} '='  'flexp.' nms{j} ';'])
end

% Parameters (setp), from char to num with the same name
nms = fieldnames(setp);
for j = 1:length(nms)
    eval([nms{j} '='  'setp.' nms{j} ';'])
end

% Steady state computation
ass       = 1;                                            
tss        = 0;                                             
chiss    = 1;
betvss  = 1;
mss     = bet;   
rss       = 1/bet-1;                                    
piiss    = 1;                                           
lamss   = 1/eta;                                                          
yss       = ass*nss^(1-alp);              
% yss       = ass*nss;                                       
css       = yss - tss;                                                     
wss      = (1-lamss)*(1-alp)*ass*nss^(-alp);   
% wss      = (1-lamss)*ass;                         
pss       = wss*css^(-sigmc)*(1-nss)^sigmn;                  
setp.psii  = pss;                                                           

% Vector of ss values. Vector must be consistent with YY and XX vectors in model.m
% y c n w m lam pii
% a t chi betv r
yyss      = [yss css nss wss mss lamss piiss];  % control variables
xxss      = [ass tss chiss betvss rss];                         % state variables
ss         = [yyss xxss];                                 % steady state values

end


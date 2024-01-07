function [ss, setp, flexp] = model_Nonsep_ss(flexp,setp)

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
G = (1-lamss)*(1-alp)/gga;
nss = G/(1+G);
yss       = ass*nss^(1-alp);              
css       = yss - tss; 
wss      = gga*css*(1-nss)^(-1);   
xiss     = css^(-sga)*(1-nss)^(gga*(1-sga));

% Vector of ss values. Vector must be consistent with YY and XX vectors in model.m
yyss       = [yss css nss wss mss lamss piiss xiss];  % control variables
xxss       = [ass tss chiss betvss rss];                       % state variables
ss         = [yyss xxss];                                 % steady state values


end


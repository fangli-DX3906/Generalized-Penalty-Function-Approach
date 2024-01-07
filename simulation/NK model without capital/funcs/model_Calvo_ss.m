function [ss, setp, flexp] = model_Calvo_ss(flexp,setp)

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
ass        = 1;                                           
tss        = 0;                                          
chiss      = 1;                                           
bss        = 1;
intss      = 1/bet-1;
sss        = 1;
vpss      = 1;
piiss      = 1;  
mss      = bet;
yss        = ass*nss;
css        = yss-tss;
x2ss     = bss*css^(-sigmc)*yss/(1-theta*bet);
x1ss     = x2ss*(eta-1)/eta;
mcss    = (1-theta*bet)*x1ss/(css^(-sigmc)*yss);
wss      = mcss*ass;
pss      = wss*css^(-sigmc)*(1-nss)^sigmn;
setp.psii = pss;

% Vector of ss values. Vector must be consistent with YY and XX vectors in model.m
yyss       = [yss css nss wss x1ss x2ss mcss piiss sss mss];  % control variables
xxss       = [ass tss chiss bss intss vpss];                       % state variables
ss         = [yyss xxss];                                 % steady state values

end


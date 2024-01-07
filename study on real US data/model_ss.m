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
ass          = 1;                                                       % tech in ss
tss           = 0;                                                       % tax in ss
chiss       = 1;                                                       % monetary policy shock parameter
betvss     = bet;                                                    % preference shock ss
rss           = 1/bet-1;                                              % net policy rate set by the monetary autority (central bank)
piiss        = 1;                                                       % inflation in ss
mss         = bet;                                                    % stochastic discount factor in ss
lamss       = 1/eta;                                                 % multiplier of the specific good yi demand
kss          = nss*((rss+del)/alp*eta/(eta-1))^(1/(alp-1)); % physical capital in ss
yss          = ass*kss^(alp)*nss^(1-alp);                  % output in ss
iiss          = del*kss;                                             % investment in ss
css          = yss - iiss - tss;                                    % consumption in ss
wss         = (1-lamss)*(1-alp)*kss^(alp)*nss^(-alp); % wage in ss
pss          = wss*css^(-sigmc)*(1-nss)^sigmn;        % weight on utility of leisure 
setp.psii  = pss;                                                    % save pss as a parameter

% Vector of ss values. Vector must be consistent with YY and XX vectors in model.m
yyss       = [yss css nss wss iiss mss lamss piiss];  % control variables
xxss       = [ass tss chiss betvss kss rss];                % state variables
ss          = [yyss xxss];                                          % steady state values

end


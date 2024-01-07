function [fyn, fxn, fypn, fxpn] = model_Nonsep(flexp,setp)

% This code solves the NK model with capital
% Code by Marco Brianti
% January 1, 2022

%Steady State
[ss,setp,flexp] = model_Nonsep_ss(flexp,setp); % 

% Parameters (flexp and setp), from char to num with the same name
nms = fieldnames(flexp);
for j = 1:length(nms)
    eval([nms{j} '='  'flexp.' nms{j} ';'])
end

nms = fieldnames(setp);
for j = 1:length(nms)
    eval([nms{j} '='  'setp.' nms{j} ';'])
end

disp(['psii = ', num2str(psii)])

% Declare Needed Symbols - see below for explanation
syms a t chi betv r 
syms a_p t_p chi_p betv_p r_p 
syms y c n w m lam pii xi
syms y_p c_p n_p w_p m_p lam_p pii_p xi_p

% %Declare X and Y vectors
XX  = [a t chi betv r]; % vector of state variables. Note al is a at t-1, and so on
XXP = [a_p t_p chi_p betv_p r_p]; % p signifies t+1 
YY  = [y c n w m lam pii xi]; % vector of control variables
YYP = [y_p c_p n_p w_p m_p lam_p pii_p xi_p]; % p signifies t+1 

%Make index variables for future use
make_index([YY,XX])
make_index([YY,XX], 2)

% Model Equations 
f(1) = - w + gga*c*(1-n)^(-1);
f(end+1) = - m + bet*(xi_p/xi);
f(end+1) = - 1 + m*(1+r_p)/pii_p;
f(end+1) = -xi + betv*c^(-sga)*(1-n)^(gga*(1-sga));
f(end+1) = - gamp*(pii-1)*pii + 1 + m*gamp*(pii_p-1)*pii_p*y_p/y - eta*lam;
f(end+1) = - w + (1-lam)*a*(1-alp)*n^(-alp); 
f(end+1) = - y + a*n^(1-alp); 
f(end+1) = - (1+r_p) + (1+r)^rhoi*((1/bet)*(pii/1)^pp)^(1-rhoi)*chi;
f(end+1) = - y + c + t + gamp/2*(pii-1)^2*y; 
f(end+1) = - log(a_p) + rhoa*log(a);
f(end+1) = - t_p + rhot*t;
f(end+1) = - log(chi_p) + rhoc*log(chi);
f(end+1) = - log(betv_p) + rhob*log(betv);

% Check Computation of Steady-State Numerically
fnum = double(subs(f, [YY XX YYP XXP], [ss, ss]));
disp('Checking steady-state equations:')
disp(fnum);

%Log-linear approx
log_var = [];
f       = subs(f, log_var, exp(log_var)); 
   
%Differentiate
fx  = jacobian(f,XX);
fy  = jacobian(f,YY);
fxp = jacobian(f,XXP); 
fyp = jacobian(f,YYP);

%Plug back into levels
fx =  subs(fx , log_var, log(log_var));
fy =  subs(fy , log_var, log(log_var));
fxp = subs(fxp, log_var, log(log_var));
fyp = subs(fyp, log_var, log(log_var));

%Compute numerical values
fxn =  double(subs(fx , [YY XX YYP XXP], [ss, ss]));
fyn =  double(subs(fy , [YY XX YYP XXP], [ss, ss]));
fxpn = double(subs(fxp, [YY XX YYP XXP], [ss, ss]));
fypn = double(subs(fyp, [YY XX YYP XXP], [ss, ss]));
end

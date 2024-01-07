function [fyn, fxn, fypn, fxpn] = model_Calvo(flexp,setp)

% This code solves the NK model with capital
% Code by Marco Brianti
% January 1, 2022

%Steady State
[ss,setp,flexp] = model_Calvo_ss(flexp,setp); % 

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
syms a t chi b int vp
syms a_p t_p chi_p b_p int_p vp_p
syms y c n w x1 x2 mc pii s m
syms y_p c_p n_p w_p x1_p x2_p mc_p pii_p s_p m_p

% %Declare X and Y vectors
XX  = [a t chi b int vp]; %betv % vector of state variables. Note al is a at t-1, and so on
XXP = [ a_p t_p chi_p b_p int_p vp_p]; % betv_p % p signifies t+1 
YY  = [y c n w x1 x2 mc pii s m]; % vector of control variables
YYP = [y_p c_p n_p w_p x1_p x2_p mc_p pii_p s_p m_p]; % p signifies t+1 

%Make index variables for future use
make_index([YY,XX])
make_index([YY,XX], 2)

% Model Equations 
f(1)        = - w + c^sigmc*psii*(1-n)^(-sigmn); 
f(end+1) = - 1 + m *(1+int_p/pii_p); 
f(end+1) = - y + c + t;
f(end+1) = - y + a*n/vp;
f(end+1) = - mc + w/a;
f(end+1) = - m + bet*(c_p/c)^(-sigmc)*(b_p/b);
f(end+1) = -1 + (1-theta)*s^(1-eta) + theta*pii^(eta-1);
f(end+1) = -vp_p + (1-theta)*s^(-eta)+theta*pii^(eta)*vp;
f(end+1) = -s+(eta/(eta-1))*x1/x2;
f(end+1) = -x1+c^(-sigmc)*mc*y*b+theta*bet*x1_p*pii_p^(eta);
f(end+1) = -x2+c^(-sigmc)*y*b+theta*bet*x2_p*pii_p^(eta-1);
f(end+1) = - log(a_p) + rhoa*log(a);    
f(end+1) = - t_p + rhot*t;              
f(end+1) = - log(chi_p) + rhoc*log(chi);                    
f(end+1) = - (1+int_p) + (1+int)^rhoi*((1/bet)*(pii/1)^pp)^(1-rhoi)*chi;                                                    
f(end+1) = -log(b_p) + rhob*log(b);

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


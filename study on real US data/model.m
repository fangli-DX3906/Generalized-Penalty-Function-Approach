function [fyn, fxn, fypn, fxpn, ss] = model(flexp,setp)

% This code solves the NK model with capital
% Code by Marco Brianti
% January 1, 2022

%Steady State
[ss,setp,flexp] = model_ss(flexp,setp); % 

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
syms a t chi betv k r
syms a_p t_p chi_p betv_p k_p r_p
syms y c n w ii m lam pii 
syms y_p c_p n_p w_p ii_p m_p lam_p pii_p

% %Declare X and Y vectors
XX  = [a t chi betv k r]; % vector of state variables. Note al is a at t-1, and so on
XXP = [a_p t_p chi_p betv_p k_p r_p]; % p signifies t+1 
YY  = [y c n w ii m lam pii]; % vector of control variables
YYP = [y_p c_p n_p w_p ii_p m_p lam_p pii_p]; % p signifies t+1 

%Make index variables for future use
make_index([YY,XX])
make_index([YY,XX], 2)

% Model Equations 
f(1)        = - w + c^sigmc*psii*(1-n)^(-sigmn);                                                    % MRS
f(end+1) = - 1 + betv*(c_p/c)^(-sigmc)*(1+r_p)/pii_p;                                         % Euler
f(end+1) = - y + a*k^(alp)*n^(1-alp);                                                                  % production function
f(end+1) = - log(a_p) + rhoa*log(a);                                                                    % TFP LOM
f(end+1) = - gamp*(pii-1)*pii + 1 + m*gamp*(pii_p-1)*pii_p*y_p/y - eta*lam;     % firm's problem FOC (p_{it})
f(end+1) = - 1 + m*(1-del+(1-lam_p)*a_p*alp*k_p^(alp-1)*n_p^(1-alp)-0.5*alp*gamp*(pii_p-1)^2*a_p*k_p^(alp-1)*n_p^(1-alp));             
% firm's problem FOC (k_{it+1})
f(end+1) = - w + (1-lam)*(1-alp)*a*k^(alp)*n^(-alp)-0.5*gamp*(pii-1)^2*(1-alp)*a*k^alp*n^(-alp);
% FOC for wage, labor demand
f(end+1) = - (1+r_p) + (1+r)^rhoi*((1/bet)*(pii/1)^pp)^(1-rhoi)*chi;                       % firm's problem FOC (b_{it}) + monetray policy
f(end+1) = - log(chi_p) + rhoc*log(chi);                                                               % monetary policy LOM
f(end+1) = - ii + k_p - (1-del)*k;                                                                          % capital LOM
f(end+1) = - y + c + ii + t + 0.5*gamp*(pii-1)^2*y;                                               % market clear condition 
f(end+1) = - t_p + rhot*t;                                                                                     % tax LOM
f(end+1) = - m + betv*(c_p/c)^(-sigmc);                                                               % definition of SDF
f(end+1) = -log(betv_p/bet) + rhob*log(betv/bet);                                                  % preference shock LOM


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

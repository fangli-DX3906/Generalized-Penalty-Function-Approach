function [delta_opt, GAM] = solver_uniform(pemodel, vars, horiz, sign, sShockPos, isCC, precision)

if isempty(precision)
    precision = 10^(-5);  
end

options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

disp('Iterations Starts...')
disp('***********************************************')

objFunc = @(delta) objPFA_unifrom(pemodel, vars, horiz, sign, sShockPos, isCC, delta);
[delta_opt] = fmincon(objFunc, [0;0], [], [], [], [], [], [], @(delta) constraint_uniform(pemodel, vars, horiz, sign, sShockPos, isCC, delta), options);

GAM = identifyPFA_deltas(pemodel, vars, horiz, sign, delta_opt, sShockPos, isCC);

end


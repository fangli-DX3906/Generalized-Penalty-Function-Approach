function solution = solver_uniform_theta(pemodel, vars, horiz, sign, sShockPos, isCC, precision, theta)

if isempty(precision)
    precision = 10^(-6);  
end

options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

disp('Iterations Starts...')
disp('***********************************************')


objFunc = @(delta) objPFA_unifrom_theta(pemodel, vars, horiz, sign, sShockPos, isCC, delta, theta);
[delta_opt] = fmincon(objFunc, [0;0], [], [], [], [], [], [], @(delta) constraint_uniform(pemodel, vars, horiz, sign, sShockPos, isCC, delta), options);
solution.GAM = identifyPFA_deltas(pemodel, vars, horiz, sign, delta_opt, sShockPos, isCC);
solution.theta = theta;
solution.delta = delta_opt;


end


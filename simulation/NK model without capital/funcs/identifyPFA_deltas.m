function GAM = identifyPFA_deltas(pemodel, vars, horiz, sign, deltas, sShockPos, isCC)
% this function allows different deltas for different PFs 

model = pemodel.basicModel;
nvar = model.NumSeries;
I_n = eye(nvar);

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

GAM = [];
for i=1:size(sShockPos, 2)
    
    objVar = vars{i};
    objHoriz = horiz(i,:);
    objSign = sign(i,:);
    delta = deltas(i);
    
    objFunc = @(gamma) objPFA(pemodel, objVar, objHoriz, objSign, gamma, delta, isCC);
    [gamma_opt] = fmincon(objFunc, I_n(:, sShockPos(i)), [], [], [], [], [], [], ...
                                       @(gamma) constraint_orthogonality(gamma), options);
                                   
    if (gamma_opt'*gamma_opt - 1)^2 > 10^(-10) 
      warning('The problem is not consistent with the constraints.')
    end   
    
    GAM = [GAM gamma_opt];
end

end


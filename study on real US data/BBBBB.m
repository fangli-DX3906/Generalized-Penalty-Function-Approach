function [GAMMA, GAM] = BBBBB(pemodel, A,  vars, horiz, sign, delta, sShockPos, isCC)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Suppose we have S set of infomation, and need S steops for identification
% objVAR should be 1 x S vector contains the variables to max SRhorizon impact in each step
% srHoriz should be 1 x S vector contains short run horizon of the maximization in each step
% initIV should be S x X vector with each row been IV used in each step
% initIVSign should be S x X matrix containing the sign informtion
% srHorizIV is the solo number for IV
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  We could make a struct that stores the information !
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

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
    
    objFunc = @(gamma) AAAAA(pemodel, A, objVar, objHoriz, objSign, gamma, delta, isCC);
    [gamma_opt] = fmincon(objFunc, I_n(:, sShockPos(i)), [], [], [], [], [], [], ...
                                       @(gamma) constraint_orthogonality(gamma), options);
                                   
    if (gamma_opt'*gamma_opt - 1)^2 > 10^(-10) 
      warning('The problem is not consistent with the constraints.')
    end   
    
    GAM = [GAM gamma_opt];
    
end

GAMMA = [GAM null(GAM')]; 

end


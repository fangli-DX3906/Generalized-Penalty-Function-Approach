function [GAMMA, delta_final, model] = solver(pemodel, vars, horiz, sign, delta, sShockPos, isCC, precision, mkInfo)

if isempty(precision)
    precision = 10^(-8);   % 10^(-7)
end

objgam = 1;
threshold = 0.1;      

% using DEN case
% add = 0.00048;

% not using DEN case
% add = 0.05;   faster!
add = 0.01;

disp('Iterations Starts...')
disp('***********************************************')
outerJ = 0;
innerJ = 0;

while objgam >= precision
%     outerJ = outerJ + 1;
%     outerJ
      while objgam >= threshold
%             innerJ = innerJ + 1;
%             innerJ
            warning on
            [GAMMA, GAM] = identifyPFA(pemodel, vars, horiz, sign, delta, sShockPos, isCC);
            objgam = GAM(:, sShockPos(1))'*GAM(:, sShockPos(2));
            delta  = delta + add;
      end
      threshold = threshold/10;
      add = add/10;
      disp(['In this epoch, the target inner product is ', num2str(objgam)])
      disp(['In this epoch, delta is ', num2str(delta)])
      disp('***********************************************')
end

if abs(objgam) == objgam
      fprintf('\n')
      disp('Identification achieved!')   
      disp(['In this epoch, delta is ', num2str(delta)])
      disp(['In this epoch, the target inner product is ', num2str(objgam)])
      fprintf('\n')
else
      fprintf('\n')
      warning(['The target inner product is ' , num2str(objgam) ])
      fprintf('\n')
end

if mkInfo == 1
    pemodel.auxMats.PFMats = pemodel.auxMats.cholVC * GAMMA;
    info.var = vars;
    info.horiz = horiz;
    info.sign = sign;
    info.shockpos = sShockPos;
    info.isCC = isCC;
    pemodel.PFidentifyInfo = info;
end

model = pemodel;

% the following lines are added
delta_final = delta;

end


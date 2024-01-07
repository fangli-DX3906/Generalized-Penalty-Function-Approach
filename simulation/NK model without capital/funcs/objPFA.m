function objPFA = objPFA(pemodel, vars, horiz, sign, gamma, delta, isCC)
% objVAR is the variable to max SRhorizon impact
% initIV is the variable(s) we impose penalty for sign
% SRhorizon is the short run horizon of the maximization of the impact of objective variable
% gamma: is a (nvar x 1) matrix from the orthogonal matrix

whos = pemodel.basicModel.SeriesNames;
A = pemodel.auxMats.cholVC;

if isCC
    deltaLi = [1-delta, delta];
else
    deltaLi = [1, delta];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delta = 0;               
% add = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std = diag(pemodel.auxMats.stdepsvec);

pos = [];
for i=1:size(vars,2)
    pos = [pos find(whos==vars{i})];
end

penalty = [];
for i=1:size(pos,2)
    IRFs_temp = calcOrthoIRFs(pemodel, A, horiz(i), gamma);
    penalty = [penalty sum(IRFs_temp(pos(i), :))];
end

scalarFEV = std(pos , pos);
% objPFA = penalty .* deltaLi * inv(scalarFEV) * sign';
objPFA = penalty .* deltaLi * sign';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delta = 0.988;
% add = 0.00048;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STD = pemodel.auxMats.stdepsvec;
% 
% pos = [];
% for i=1:size(vars,2)
%     pos = [pos find(whos==vars{i})];
% end
% 
% IRFs_temp1 = calcOrthoIRFs(pemodel, A, horiz(1), gamma);
% IRFs_temp2 = calcOrthoIRFs(pemodel, A, horiz(2), gamma);
% xxxxx = calcOrthoIRFs(pemodel, A, horiz(1), NaN);
% DEN = sum(sum(xxxxx(pos(1),:,:).^2));
% penalty = [sum(IRFs_temp1(pos(1), :))/DEN  sum(IRFs_temp2(pos(2), :))/STD(pos(2))];
% objPFA = penalty .* deltaLi * sign';

objPFA = -objPFA;

end


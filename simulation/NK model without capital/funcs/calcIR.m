function IRFs = calcIR(pemodel, resData, irType, irHoriz, isCI, ciMethod, alpha)
% irMethod: the method of calculating the IRF  
%                 'general' --- general method (built in method)
%                 'gpfa' --- general panelty function approach (Brianti, 2021)
%  irType: only works when irMethod == 'general', 'generalized' or 'orthogonalized'
%  IRFs: a struct stores the necessary info and pass to plotIR for plotting IRs 

rng('default');
rng(12);

model = pemodel.basicModel;
    
if isCI == 1
    if ciMethod == 'mcs'   % MCS 
        [IR, lowerCI, upperCI] = irf(model, 'Method', irType, 'NumObs', irHoriz, 'NumPaths', 200, 'SampleSize', 100, 'Confidence', alpha);
    else   % bootstrapping
        [IR, lowerCI, upperCI] = irf(model, 'Method', irType, 'NumObs', irHoriz, 'E', resData, 'NumPaths', 100, 'SampleSize', 100, 'Confidence', alpha);
    end
    info.alpha = alpha; 
else
    IR = irf(model, 'Method', irType, 'NumObs', irHoriz);
    upperCI = NaN;
    lowerCI = NaN;
end
    
info.h = irHoriz;
info.whos = model.SeriesNames;
info.model = model.Description;
IRFs.irf = IR;
IRFs.uci = upperCI;
IRFs.lci = lowerCI;
IRFs.info = info;

end
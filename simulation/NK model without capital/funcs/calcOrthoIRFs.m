function IRFs = calcOrthoIRFs(peModel, A, irHoriz, initInput)
%  This is a lightweight function calculates the IRF (without CIs)
%  initInput = NaN, use cols in an identity matrix for the initial shock
%  initInput = same dimension vector, use this vector as initial shock
% if mkStruct ==1 then IRFs will be a struct

model = peModel.basicModel;
nlags = model.P;
nvar = model.NumSeries;
nshocks = model.NumSeries;

% if isPF == 1
%     A = peModel.auxMats.PFMats;
% else
%     A = peModel.auxMats.cholVC;
% end

B = peModel.auxMats.Bpl_ev(1:end-1,:);
I_n = eye(nvar);


if max(size(initInput)) == nvar  % case of providing the initial shock
    
        IRFs = zeros(nvar, irHoriz);
        
        shocks = initInput;
        IRFs(:,1) = A*shocks;
        F = [IRFs(:,1)' zeros(1,(nlags-1)*nvar)];
        
        for k = 2:irHoriz
            IRFs(:,k) = F*B;
            F = [IRFs(:,k)' F(1:end-nvar)];
        end
        
else  % case of not providing the initial shock
    
    IRFs = zeros(nvar, irHoriz, nvar);
   
    for i_shock = 1:nshocks
        
        shocks = I_n(:, i_shock);
        IRFs(:, 1, i_shock) = A * shocks;
        F = [IRFs(:, 1, i_shock)' zeros(1,(nlags-1)*nvar)];
        
        for k = 2:irHoriz
            IRFs(:, k, i_shock) = F * B;
            F = [IRFs(:,k,i_shock)' F(1:end-nvar)];
        end
        
    end
end

end


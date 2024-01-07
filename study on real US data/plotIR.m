function plotIR(pemodel, IRFs, selectVars, selectShocks, isGrid, shockNames)
% IRF should be an instance

if ~isstruct(IRFs)
    IRFs = mkIRFStruct(pemodel, IRFs, [], [], []);
end

varList = [];
% to find the endog var's positions 
for i=1:size(selectVars,2)
    varList = [varList find(IRFs.info.whos==selectVars(i))];
end

plotcount = 1;
nShocks = size(selectShocks, 2);
nVars = size(varList, 2);
plot = sum(sum(sum(~isnan(IRFs.uci))));

% for vv=1:size(varList,2)
%     
%     yaxmin = 1.1*min(min(IRFs.irf(:, :, varList(vv))));
%     yaxmax = 1.1*max(max(IRFs.irf(:, :, varList(vv))));
%     for ss=1:nShocks
%         
%         subplot(nVars, nShocks, plotcount);
%         IR = IRFs.irf(:, selectShocks(ss), varList(vv));
%         if plot
%             CI = [IRFs.lci(:, selectShocks(ss), varList(vv)) IRFs.uci(:, selectShocks(ss), varList(vv))];
%         else
%             CI = NaN;
%         end
%         % invoke plotOneIR        
%         plotOneIR(IR, CI, IRFs.info.h, selectVars(vv), shockNames{ss}, isGrid);
% %         axis( [0 (IRFs.info.h) yaxmin yaxmax] );   
%         axis tight;
%         plotcount = plotcount+1;  
%         
%     end
% end

%%%%%%%%%
for ss=1:nShocks
    
    for vv=1:size(varList,2)
        
        subplot(nShocks, nVars, plotcount);
        IR = IRFs.irf(:, selectShocks(ss), varList(vv));
        
        if plot
            CI = [IRFs.lci(:, selectShocks(ss), varList(vv)) IRFs.uci(:, selectShocks(ss), varList(vv))];
        else
            CI = NaN;
        end
        % invoke plotOneIR
        if ss==1
            plotOneIR(IR, CI, IRFs.info.h, selectVars(vv), shockNames{ss}, isGrid);
        else
            plotOneIR(IR, CI, IRFs.info.h, [], [], isGrid);
        end
%         axis( [0 (IRFs.info.h) yaxmin yaxmax] );   
        axis tight;
        plotcount = plotcount+1;  
        
    end
end
%%%%%%%%%

if plot
    legend({"IRF" [num2str(100*IRFs.info.alpha) '% Confidence Interval']}, 'Location', 'south')
% else
%     legend("IRF")
end

end


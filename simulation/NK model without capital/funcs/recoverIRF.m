function IRF_new = recoverIRF(nShocks, nVar, nAux, IRFr, As, plist, h, auxcon)

IRF_new = [];

for s=1:nShocks 
    % pick out demand and supply entries
    if s==1
        sv = 1:nVar;              % demand shock
    else
        sv = nVar+1:2*nVar;  % supply shock
    end
    
    for j = 1:nAux
        Aa = As{j};
        par = plist(j);
        
        for hh = 1: h+1
            
            if hh <= par+1
                aaa =  IRFr(sv, 1:hh);
                bbb = [];
                for idx=1:hh
                    bbb = [aaa(:, idx); bbb];
                end
                if auxcon
                    irff = Aa(:, 1:1+nVar*hh) * [1; bbb];
                else
                    irff = Aa(:, 1:nVar*hh) * bbb;
                end
                IRF_new = [IRF_new, irff];
            else
                aaa =  IRFr(sv, hh-par:hh);
                bbb = [];
                for idx=1:par+1
                    bbb = [aaa(:, idx); bbb];
                end
                if auxcon
                    irff = Aa * [1; bbb];
                else
                    irff = Aa * bbb;
                end
                IRF_new = [IRF_new, irff];
            end  
            
        end  % end of this period   
        
    end  % end of this aux variable 
    
end % end of this shock

end


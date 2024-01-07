function optim_solutions = solver_deltas_shrinks(pemodel, vars, horiz, sign, sShockPos, isCC, precision)

if isempty(precision)
    precision = 10^(-5);   
end

objgam_init = 1;
threshold_init = 0.1; 
add_init = 0.05;

disp('Iterations Starts...')
disp('*********************************************************************')
delta1_list = 0:0.01:1;
optim_solutions = {};

for d = 1:length(delta1_list)
    threshold = threshold_init;
    add = add_init;
    objgam = objgam_init;
    delta1 = delta1_list(d);
    delta2 = 0;
    flag =1;
    disp(['A new epoch, delta 1= ', num2str(delta1)])
    while objgam >= precision
        disp('   another batch...')
        
        while objgam >=threshold
            delt = [delta1, delta2];
            delt
            GAM = identifyPFA_deltas(pemodel, vars, horiz, sign, delt, sShockPos, isCC);
            objgam = abs(GAM(:, sShockPos(1))'*GAM(:, sShockPos(2)));
            delta2  = delta2 + add;
            if delta2>1
                break
            end
        end
        
        threshold = threshold/10;
        add = add/10;
        objgam
        
        if delta2>1
            disp(['   Fail to achieve the identification in this epoch'])
            flag = 0;
            break
        end
        
    end
    
    if flag == 1
        disp(['In this epoch, delta_1 = ', num2str(delt(1)), ', delta_2 = ', num2str(delt(2))])
        disp(['In this epoch, the target inner product is ', num2str(objgam)])
        disp('***********************************************')
        solut.delta = delt;
        solut.gamma = [GAM null(GAM')];
        optim_solutions{end+1} = solut;
    end
    
end

end


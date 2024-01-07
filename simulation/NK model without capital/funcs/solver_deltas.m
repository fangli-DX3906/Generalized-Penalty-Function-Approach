function optim_solutions = solver_deltas(pemodel, vars, horiz, sign, sShockPos, isCC, precision)

if isempty(precision)
    precision = 10^(-5);  
end

% step = 0.002;
% delta1_list = 0:step:1;
% delta2_list = 0:step:1;
delta1_list =0:0.001:1;
cv = pemodel.basicModel.Covariance;

disp('Iterations Starts...')
disp('***********************************************')

optim_solutions = {};

for d1 = 1:length(delta1_list)
    delta1 = delta1_list(d1);
%     delta1= sqrt(cv(2,2))/(sqrt(cv(1,1))+sqrt(cv(2,2)));
    theta1 = 1/delta1-1;
    theta2 = (cv(1,1)+cv(1,2)*theta1)/(cv(2,1)+cv(2,2)*theta1);
    delta2 = 1/(theta2+1);
    if delta2<1 && delta2>0
        delta = [1-delta1, 1-delta2];
        GAM = identifyPFA_deltas(pemodel, vars, horiz, sign, delta, sShockPos, isCC);
        objgam = GAM(:, sShockPos(1))'*GAM(:, sShockPos(2));
        if abs(objgam)<=precision
            disp(['The target inner product is ', num2str(objgam)])
            disp(['Delta are ', num2str(1-delta1),'     ', num2str(1-delta2)])
            solut.delta = delta;
            solut.gamma = [GAM null(GAM')];
            optim_solutions{end+1} = solut;
        end
    end
end

end


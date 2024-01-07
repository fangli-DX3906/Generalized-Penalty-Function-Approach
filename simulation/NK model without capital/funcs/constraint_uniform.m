  function [c,ceq] = constraint_uniform(pemodel, vars, horiz, sign, sShockPos, isCC, delta)
       
         gamma = identifyPFA_deltas(pemodel, vars, horiz, sign, delta, sShockPos, isCC);
         
         c = [delta(1)-1;
                -delta(1);
                delta(2)-1;
                -delta(2)];
         ceq = gamma'*gamma-eye(size(gamma,2));
         
   end
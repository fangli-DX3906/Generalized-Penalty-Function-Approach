function Jacobian = calcJacobianStructuralwrtReduced(reduced_var_param, covariance_mat, delta_1, precision)

% structural param (A_0, A_+): n^2+n*(n*p+1) = n^2(p+1)+n
% reduced param (B, sigma, delta_1): n^2+n*(n*p+1)+1 = n^2(p+1)+n+1
% so jacobian is of (n^2(p+1)+n+1) x (n^2(p+1)+n) dimension.

B = reduced_var_param;
sigma = covariance_mat;
n = size(sigma, 1);
p = (size(B, 1)-1)/n;

Jacobian = zeros(n^2+n*(n*p+1),  n*(n*p+1)+n^2+1);

% (B, Sigma, delta1) ---> (Aplus, A0)
for i=1:n*(n*p+1)+n^2+1
    if i<=n*(n*p+1)
        B_vec_ = vec(B);  
        B_vec_(i,:) = B_vec_(i, :) + precision;
        B_ = reshape(B_vec_, size(B,1), []);
        [A_0, A_plus] = calcStructuralVARParams(B, sigma, delta_1);
        [A_0_, A_plus_] = calcStructuralVARParams(B_, sigma, delta_1);
    elseif i>n*(n*p+1) && i<=n*(n*p+1)+n^2
        sigma_vec_ = vec(sigma);
        if i-n*(n*p+1)==2
            sigma_vec_(i-n*(n*p+1),:) = sigma_vec_(i-n*(n*p+1), :) + precision; 
            sigma_vec_(i-n*(n*p+1)+1,:) = sigma_vec_(i-n*(n*p+1)+1, :) + precision; 
        elseif i-n*(n*p+1)==3
            sigma_vec_(i-n*(n*p+1),:) = sigma_vec_(i-n*(n*p+1), :) + precision; 
            sigma_vec_(i-n*(n*p+1)-1,:) = sigma_vec_(i-n*(n*p+1)-1, :) + precision; 
        else
            sigma_vec_(i-n*(n*p+1),:) = sigma_vec_(i-n*(n*p+1), :) + precision; 
        end
        sigma_ = reshape(sigma_vec_, size(sigma,1), []);
        [A_0, A_plus] = calcStructuralVARParams(B, sigma, delta_1);
        [A_0_, A_plus_] = calcStructuralVARParams(B, sigma_, delta_1);
    else
        [A_0, A_plus] = calcStructuralVARParams(B, sigma, delta_1);
        [A_0_, A_plus_] = calcStructuralVARParams(B, sigma, delta_1 + precision);
    end
    A = [vec(A_0); vec(A_plus)];
    A_ = [vec(A_0_); vec(A_plus_)];
    Jacobian(:, i) = (A_ - A)/precision;
end

% (Sigma, B, delta1) ---> (A0, Aplus)
% for i=1:n*(n*p+1)+n^2+1
%     if i<= n^2
%         sigma_vec_ = vec(sigma);
%         if i==2
%             sigma_vec_(i,:) = sigma_vec_(i, :) + precision; 
%             sigma_vec_(i+1,:) = sigma_vec_(i+1, :) + precision; 
%         elseif i==3
%             sigma_vec_(i,:) = sigma_vec_(i, :) + precision; 
%             sigma_vec_(i-1,:) = sigma_vec_(i-1, :) + precision; 
%         else
%             sigma_vec_(i,:) = sigma_vec_(i, :) + precision; 
%         end
%         sigma_ = reshape(sigma_vec_, size(sigma,1), []);
%         [A_0, A_plus] = calcStructuralVARParams(B, sigma, delta_1);
%         [A_0_, A_plus_] = calcStructuralVARParams(B, sigma_, delta_1);
%     elseif i>n^2 && i<=n*(n*p+1)+n^2
%         B_vec_ = vec(B);  
%         B_vec_(i-n^2,:) = B_vec_(i-n^2, :) + precision;
%         B_ = reshape(B_vec_, size(B,1), []);
%         [A_0, A_plus] = calcStructuralVARParams(B, sigma, delta_1);
%         [A_0_, A_plus_] = calcStructuralVARParams(B_, sigma, delta_1);
%     else
%         [A_0, A_plus] = calcStructuralVARParams(B, sigma, delta_1);
%         [A_0_, A_plus_] = calcStructuralVARParams(B, sigma, delta_1 + precision);
%     end
%     A = [vec(A_0); vec(A_plus)];
%     A_ = [vec(A_0_); vec(A_plus_)];
%     Jacobian(:, i) = (A_ - A)/precision;
% end

% try another: is this reasonable?
% for i=1:n*(n*p+1)+n^2
%     if i<=n*(n*p+1)
%         B_vec_ = vec(B);  
%         B_vec_(i,:) = B_vec_(i, :) + precision;
%         B_ = reshape(B_vec_, size(B,1), []);
%         [A_0, A_plus] = calcStructuralVARParams(B, sigma, delta_1);
%         [A_0_, A_plus_] = calcStructuralVARParams(B_, sigma, delta_1);
%     elseif i>n*(n*p+1) && i<=n*(n*p+1)+n^2-1
%         sigma_vec_ = vec(sigma);
%         if i-n*(n*p+1)==2
%             sigma_vec_(i-n*(n*p+1),:) = sigma_vec_(i-n*(n*p+1), :) + precision; 
%             sigma_vec_(i-n*(n*p+1)+1,:) = sigma_vec_(i-n*(n*p+1)+1, :) + precision; 
%         else
%             sigma_vec_(i-n*(n*p+1),:) = sigma_vec_(i-n*(n*p+1), :) + precision; 
%         end
%         sigma_ = reshape(sigma_vec_, size(sigma,1), []);
%         [A_0, A_plus] = calcStructuralVARParams(B, sigma, delta_1);
%         [A_0_, A_plus_] = calcStructuralVARParams(B, sigma_, delta_1);
%     else
%         [A_0, A_plus] = calcStructuralVARParams(B, sigma, delta_1);
%         [A_0_, A_plus_] = calcStructuralVARParams(B, sigma, delta_1 + precision);
%     end
%     A = [vec(A_0); vec(A_plus)];
%     A_ = [vec(A_0_); vec(A_plus_)];
%     Jacobian(:, i) = (A_ - A)/precision;
% end

% just to transpose it.
Jacobian = Jacobian';

end


function J3 = calcJ3(n, H, p, A_0, A_plus, covariance_mat, companion_mat)

sigma = covariance_mat;
IRF = calcReducedIRF(companion_mat, sigma, p, H);
A_plus = A_plus(2:end, :);   % remove the intercept

K = makeCommutationMat(n,n);
M = [eye(n), zeros(n, n*(p-1))]';    % size np x n
F = [];
for i=1:p
    Ai = A_plus((i-1)*n+1:i*n, :);
    Bi = Ai * inv(A_0);
    F = [F Bi'];
end

% F matrix seems to be the companion matrix defined in Kilian and Lutkepohl, 
% but in Inuoe and Kilian (2022) F is defined as [F; zeros(n*(p-1), n) eye(n*(p-1))] (eq 27).
% This is not the same companion matrix. Also, in Inuoe and Kilian (2022)
% they claim this F matrix follows Eq. (A.25) in Cheng et al. (2020) and Cheng et al says F is the companion matrix.
% So maybe this is a typo in Inuoe and Kilian (2022), here we use the correct companion matrix for now.
F = [F; eye(n*(p-1)) zeros(n*(p-1), n)];     %  size: np x np

if H<=p
    
    J3 = eye(n^2*(H+1));
    J3(1:n^2, 1:n^2) = -kron(inv(A_0), inv(A_0')) * K;
    
    for i=3:H+1
        X = calcXPhiB(M, F, IRF, K, n, p, i-1);   % size: n^2 x n^2p
        J3( (i-1)*n^2+1:i*n^2, n^2+1:(i-1)*n^2 ) = X(:, 1:(i-2)*n^2);
    end
    
else
    
    J3 = eye(n^2*(H+1), n^2*(p+1));
    J3(1:n^2, 1:n^2) = -kron(inv(A_0), inv(A_0')) * K;
    for i=3:H+1
        X = calcXPhiB(M, F, IRF, K, n, p, i-1);   % size: n^2 x n^2p
        if i<=p
           J3( (i-1)*n^2+1:i*n^2, n^2+1:(i-1)*n^2 ) = X(:, 1:(i-2)*n^2);
        else
           J3( (i-1)*n^2+1:i*n^2, n^2+1:end ) = X;
        end
    end
    
end

end
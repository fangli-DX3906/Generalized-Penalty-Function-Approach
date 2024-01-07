function J4 = calcJ4(n, H, p, A_0, A_plus)

K = makeCommutationMat(n,n);
m = min(p, H);
A_plus = A_plus(2:end, :);

J4 = zeros(n^2*(m+1));
J4(1:n^2, 1:n^2)=eye(n^2);

for i=2:m+1
    Ai = A_plus((i-2)*n+1:(i-1)*n, :);
    J4((i-1)*n^2+1: i*n^2, 1:n^2) = -kron((Ai*inv(A_0)), inv(A_0')) * K;
    J4((i-1)*n^2+1: i*n^2, (i-1)*n^2+1: i*n^2) = kron(eye(n), inv(A_0')) * K;   
end

end


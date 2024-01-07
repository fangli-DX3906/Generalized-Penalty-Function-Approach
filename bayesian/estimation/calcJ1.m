function J1 = calcJ1(cumulative_vars, n, H)

if isempty(cumulative_vars)
    J1 = eye(n^2*(H+1));
else
    S = zeros(n,1);
    S(cumulative_vars) = 1;
    S = diag(S);
    J1 = eye(n^2*(H+1));
    for i=2:H+1
        s = (i-1)*n^2+1;
        e = i*n^2;
        for j=1:i-1
            J1(s:e, (j-1)*n^2+1:j*n^2) = kron(eye(n), S);
        end
    end
end

end


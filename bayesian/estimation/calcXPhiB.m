function X = calcXPhiB(M, F, IRF, K, n, p, i)

X = 0;
for k=0:i-1
    phi = IRF(:, i-1-k+1);
    X = X + kron(M'*(F^k)', reshape(phi, n, [])) * kron(eye(p), K);
end

end


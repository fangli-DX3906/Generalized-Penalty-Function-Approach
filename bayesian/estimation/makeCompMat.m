function comp_mat = makeCompMat(reduced_var_param, n, p)

B = reduced_var_param;
B = B(2:end, :);
comp_mat = [B'; eye(n*p-n, n*p)];

end


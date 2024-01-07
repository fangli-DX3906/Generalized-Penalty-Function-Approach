function Q = calcGPFARotation(covariance_mat, delta_1)

e1 = [1, 0];
e2 = [0, 1];
sigma = covariance_mat;

if delta_1==0
    theta_1 = eps;
else
    theta_1 = 1/delta_1-1;
end

theta_2 = (sigma(1,2) + sigma(2,2)*theta_1) / (sigma(1,1) + sigma(1,2)*theta_1);

q = zeros(2,2);

% according to GPFA_algorithm_3.pdf  page4
q(1,1) = 1 / sqrt( (e1+theta_1*e2) * sigma * (e1'+theta_1*e2') );
q(1,2) = theta_2 / sqrt( (theta_2*e1-e2) * sigma * (theta_2*e1'-e2') );
q(2,1) = theta_1 / sqrt( (e1+theta_1*e2) * sigma * (e1'+theta_1*e2') );
q(2,2) = -1 / sqrt( (theta_2*e1-e2) * sigma * (theta_2*e1'-e2') );

Q = chol(sigma) * q;

end


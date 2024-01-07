function Jacobian = calcJacobianIRFwrtStructural(J1, J2, J3, J4)

if size(J1, 2)~=size(J2, 1) || size(J2, 2)~=size(J3, 1) || size(J3, 2)~=size(J4, 1)
    disp('dimension imcompatible!');
    Jacobian = nan;
else
    Jacobian = J1*J2*J3*J4;
end

end


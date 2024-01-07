function optimOrder = auxorderSelect(sdata, ydata, maxOrder, auxcon)

nVar = size(sdata, 2);
nNum = size(sdata, 1);
nNVar = size(ydata, 2);
aicList = [];

for i=1:maxOrder
     [~, SIGMAaux]  = auxreg(sdata, ydata, i, auxcon);
    VC_eps = SIGMAaux(1:nVar, 1:nVar);
    aicList = [aicList, log(det(VC_eps)) + (2*nVar*nNVar*i) / nNum];     % no intercept
end

[~, optimOrder] = min(aicList);

end


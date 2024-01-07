function optimOrder = auxorderSelect(ydata, sdata, maxOrder)

nVar = size(sdata, 2);
nNum = size(sdata, 1);
aicList = [];


for i=1:maxOrder
     [~, SIGMAaux]  = olsvaraux(sdata, ydata, i);
    VC_eps = SIGMAaux(1:nVar, 1:nVar);
    aicList = [aicList, log(det(VC_eps)) + (2*nVar^2*i) / nNum];
end

[~, optimOrder] = min(aicList);

end


function optimOrder = orderSelect(ydata, maxOrder)

nVar = size(ydata, 2);
nNum = size(ydata, 1);
aicList = [];

logList = [];
numList = [];

for i=1:maxOrder
    [~,SIGMA,~,~,~] = olsvarc(ydata, i);
    VC_eps = SIGMA(1:nVar, 1:nVar);
    aicList = [aicList, log(det(VC_eps)) + (2*nVar^2*i) / nNum];
end

[~, optimOrder] = min(aicList);

end


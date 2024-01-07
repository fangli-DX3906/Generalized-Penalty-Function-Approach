function optimOrder = orderSelect(ydata, maxOrder)

nVar = size(ydata, 2);
nNum = size(ydata, 1);
aicList = [];

logList = [];
numList = [];

for i=1:maxOrder
    % V0
%     md = varm(nVar, i);
%     [~,~,thisLog,~] = estimate(md, ydata);
%     logList = [logList thisLog];
%     numList = [numList nVar^2*i];
%     aicList = aicbic(logList, numList);
    % V1
%     [~, ~, ~, ~, VC_eps, ~] = estim(ydata, i, 1, 0);
    % V2
    [~,SIGMA,~,~,~] = olsvarc(ydata, i);
    VC_eps = SIGMA(1:nVar, 1:nVar);
    aicList = [aicList, log(det(VC_eps)) + (2*nVar^2*i) / nNum];
end

[~, optimOrder] = min(aicList);

end


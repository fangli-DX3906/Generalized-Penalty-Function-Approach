function [Aset, orderList] = olsvaraux(s, y, iscon, maxOrder)

orderList = [];
Aset = cell(size(s,2),1);

for n = 1: size(s,2)
    p = auxorderSelect(s(:,n), y, maxOrder, iscon);
    orderList = [orderList, p];
    [A, ~] = auxreg(s(:, n), y, p, iscon);
    Aset{n} = A;
end
    
end

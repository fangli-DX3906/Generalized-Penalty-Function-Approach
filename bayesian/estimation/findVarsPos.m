function result = findVarsPos(Alist, Blist)
 
result = nan(1, size(Alist,2));

for i=1:size(Alist,2)
    for j=1:size(Blist, 2)
        if string(Alist{i})==string(Blist{j})
            result(i)=j;
        end
    end
end

if sum(isnan(result))
    disp(['Some variables are not found']);
    result = result(~isnan(result));
end

end



function [data, tp, tp2] = truncate_data(data)

n_var_system = size(data,2);
threshold    = -1/eps;
for i_var_system = 1:n_var_system
      loc(i_var_system) = find(data(:,i_var_system) > threshold, 1);
end
tp   = max(loc);
data = [data(tp:end,:); NaN(1,n_var_system)];

loc2 = zeros(1,n_var_system);
for i_var_system2 = 1:n_var_system
      loc2(1,i_var_system2) = find(isnan(data(:,i_var_system2)),1);
end
tp2  = min(loc2) - 1;
data = data(1:tp2,:);
tp2  = tp + tp2 - 1;

end
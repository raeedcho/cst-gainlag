function q = calc_tangling(x,dx)
% CALC_TANGLING calculates tanging of series x, given the derivatives dx
% Inputs:
%   x - TxP matrix, where rows are time observations and columns are features
%   dx - TxP matrix of time derivatives of x

norm_const = 1e-6;
q = zeros(size(x,1),1);
for i = 1:size(x,1)
    temp_tang = zeros(size(x,1),1);
    for j = 1:size(x,1)
        if i ~= j
            temp_tang(j) = norm(dx(i,:)-dx(j,:))^2/(norm(x(i,:)-x(j,:))^2+norm_const);
        end
    end
    temp_tang(i) = [];
    
    if 0
        q(i) = max(temp_tang);
    else
        q(i) = prctile(temp_tang,99.99);
    end
end
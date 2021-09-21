function cost = cart_pole_simpson_cost(var_list,n,h)
% Set the objective function for cart pole problem
% var_list: decision variables lists to optimize, cost: value to minimize
u_list=var_list((6*n+1):(7*n));
% Cost - Hermite-Simpson
cost=h/6*(u_list(1)^2+u_list(end)^2+4*sum(u_list(2:2:(end-1)).^2)+2*sum(u_list(3:2:(end-2)).^2));
%cost=mean(u_list.^2);
end
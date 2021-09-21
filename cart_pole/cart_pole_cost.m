function cost = cart_pole_cost(var_list,n,tf)
% Set the objective function for cart pole problem
% var_list: decision variables lists to optimize, cost: value to minimize
h=tf/(n-1);
u_list=var_list((6*n+1):(7*n));
cost=h/2*(u_list(1)^2+u_list(end)^2+2*sum(u_list(2:(end-1)).^2));
%cost=mean(u_list.^2);
end
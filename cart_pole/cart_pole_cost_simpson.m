function cost = cart_pole_simpson_timeopt_cost(var_list,n)
% Set the objective function for cart pole problem
% var_list: decision variables lists to optimize, cost: value to minimize
n_seg=(n-1)/2;
T=var_list(7*n+1);
u_list=var_list((6*n+1):(7*n));
effort_avg=1/(6*n_seg)*(u_list(1)^2+u_list(end)^2+4*sum(u_list(2:2:(end-1)).^2)+2*sum(u_list(3:2:(end-2)).^2));
weight_effort=1; 
weight_T=50;
cost=weight_effort*effort_avg + weight_T*T;
end
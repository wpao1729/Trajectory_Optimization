function cost = compute_cost(t_list,xf_list)
% Set the objective function for block-move problem
% xf_list: path to optimize, cost: value to minimize
f_list=xf_list((length(xf_list)/2+1):end);
h=(t_list(end)-t_list(1))/(length(t_list)-1);
u_list=(f_list(2:end)-f_list(1:(end-1)))/h;
cost=h*sum(u_list.^2);
end
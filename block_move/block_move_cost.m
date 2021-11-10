function cost = block_move_cost(var_list,tf,n)
% Set the objective function for block-move problem

h=tf/(n-1);
u_list=var_list((3*n+1):4*n);
effort_list=u_list.^2;
% Compute cost by integral of force squared over time, trapezoidal rule
cost=h/2*(effort_list(1)+effort_list(end)+2*sum(effort_list(2:(end-1))));

end
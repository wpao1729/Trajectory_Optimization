function cost = step_strike_cost(var_list,n,tf)

h=tf/(n-1);
u_list=var_list((55*n+1):(59*n));
u_mat=reshape(u_list,4,n);
effort_list=sum(u_mat.^2);
cost=h/2*(effort_list(1)+effort_list(end)+2*sum(effort_list(2:(end-1))));

end
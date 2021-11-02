function cost = single_stance_cost(var_list,n,tf)

h=tf/(n-1);
u_list=var_list((55*n+1):(60*n));
u_mat=reshape(u_list,5,n);
effort_list=sum(u_mat.^2);
cost=h/2*(effort_list(1)+effort_list(end)+2*sum(effort_list(2:(end-1))));

end
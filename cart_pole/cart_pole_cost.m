function cost = cart_pole_cost(var_list,n,simpson)
% Set the objective function for cart pole problem

T=var_list(7*n+1);
h=(T-0)/(n-1); % time step size
u_list=var_list((6*n+1):(7*n));

if simpson
    % Integral of force over time by simpson rule
    cost=h/3*(u_list(1)^2+u_list(end)^2+4*sum(u_list(2:2:(end-1)).^2)+2*sum(u_list(3:2:(end-2)).^2));

    % effort_avg=1/(6*n_seg)*(u_list(1)^2+u_list(end)^2+4*sum(u_list(2:2:(end-1)).^2)+2*sum(u_list(3:2:(end-2)).^2));
    % weight_effort=1; 
    % weight_T=50;
    % cost=weight_effort*effort_avg + weight_T*T;

else
    % Integral of force over time by trapezoidal rule
    cost=h/2*(u_list(1)^2+u_list(end)^2+2*sum(u_list(2:(end-1)).^2));

    % effort_avg=(u_list(1)^2+u_list(end)^2+2*sum(u_list(2:(end-1)).^2))/2/(n-1);
    % weight_effort=1; 
    % weight_T=50;
    % cost=weight_effort*effort_avg + weight_T*T;
end

end
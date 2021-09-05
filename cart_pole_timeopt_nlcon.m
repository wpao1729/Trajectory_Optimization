function [c,ceq] = cart_pole_timeopt_nlcon(var_list,n,mc,mp,L,g)
% Nonlinear Constraints for the cart pole problem w/ time optimization

x_list=var_list(1:n); v_list=var_list((n+1):(2*n)); a_list=var_list((2*n+1):(3*n)); 
theta_list=var_list((3*n+1):(4*n)); omega_list=var_list((4*n+1):(5*n)); alpha_list=var_list((5*n+1):(6*n)); 
u_list=var_list((6*n+1):(7*n)); T=var_list(7*n+1);

% System Dynamics Constraints
ceq1=u_list-(mc+mp)*a_list+mp*L*sin(theta_list).*alpha_list+mp*L*cos(theta_list).*(omega_list.^2);
ceq2=L*alpha_list-a_list.*sin(theta_list)+g*cos(theta_list);

% Collocation constraints
h=T/(n-1);
ceq3=x_list(2:end)-x_list(1:(end-1))-h/2*(v_list(1:(end-1))+v_list(2:end));
ceq4=v_list(2:end)-v_list(1:(end-1))-h/2*(a_list(1:(end-1))+a_list(2:end));
ceq5=theta_list(2:end)-theta_list(1:(end-1))-h/2*(omega_list(1:(end-1))+omega_list(2:end));
ceq6=omega_list(2:end)-omega_list(1:(end-1))-h/2*(alpha_list(1:(end-1))+alpha_list(2:end));

ceq=[ceq1,ceq2,ceq3,ceq4,ceq5,ceq6];
c=0; % No inequal constraints
end
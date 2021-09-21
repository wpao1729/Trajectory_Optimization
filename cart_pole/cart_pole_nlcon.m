function [c,ceq] = cart_pole_nlcon(var_list,n,mc,mp,L,g)
% Nonlinear Constraints for the cart pole problem - System Dynamics
a_list=var_list((2*n+1):(3*n)); theta_list=var_list((3*n+1):(4*n)); omega_list=var_list((4*n+1):(5*n)); alpha_list=var_list((5*n+1):(6*n)); u_list=var_list((6*n+1):(7*n));
ceq1=u_list-(mc+mp)*a_list+mp*L*sin(theta_list).*alpha_list+mp*L*cos(theta_list).*(omega_list.^2);
ceq2=L*alpha_list-a_list.*sin(theta_list)+g*cos(theta_list);
ceq=[ceq1,ceq2];
c=0; % No inequal constraints
end
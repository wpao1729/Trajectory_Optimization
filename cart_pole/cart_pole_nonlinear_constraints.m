function [c,ceq] = cart_pole_nonlinear_constraints(var_list,n,m_cart,m_pole,L,g,simpson)
% Nonlinear Constraints for the cart pole problem

ceq=[]; % all entries of ceq are constrained to be =0
c=[]; % all entries of c are constrained to be <=0

% Break down the decision variable list
x_list=var_list(1:n); 
v_list=var_list((n+1):(2*n)); 
a_list=var_list((2*n+1):(3*n)); 
theta_list=var_list((3*n+1):(4*n)); 
omega_list=var_list((4*n+1):(5*n)); 
alpha_list=var_list((5*n+1):(6*n)); 
u_list=var_list((6*n+1):(7*n)); 
T=var_list(7*n+1);

% Kinematics and Dynamics
ceq1=u_list-(m_cart+m_pole)*a_list+m_pole*L*sin(theta_list).*alpha_list+m_pole*L*cos(theta_list).*(omega_list.^2);
ceq2=L*alpha_list-a_list.*sin(theta_list)+g*cos(theta_list);

% Collocation constraints (nonlinear due to h dependent on T)
if simpson
    % Hermite-simpson collocation
    n_seg=(n-1)/2;
    h=T/n_seg;
    ceq3=zeros(1,n_seg); ceq4=zeros(1,n_seg); ceq5=zeros(1,n_seg); ceq6=zeros(1,n_seg); ceq7=zeros(1,n_seg); ceq8=zeros(1,n_seg); ceq9=zeros(1,n_seg); ceq10=zeros(1,n_seg);
    for i=1:n_seg
        ceq3(i)=x_list(2*i+1)-x_list(2*i-1)-h/6*(v_list(2*i-1)+4*v_list(2*i)+v_list(2*i+1));
        ceq4(i)=x_list(2*i)-1/2*(x_list(2*i-1)+x_list(2*i+1))-h/8*(v_list(2*i-1)-v_list(2*i+1));
        ceq5(i)=v_list(2*i+1)-v_list(2*i-1)-h/6*(a_list(2*i-1)+4*a_list(2*i)+a_list(2*i+1));
        ceq6(i)=v_list(2*i)-1/2*(v_list(2*i-1)+v_list(2*i+1))-h/8*(a_list(2*i-1)-a_list(2*i+1));
        ceq7(i)=theta_list(2*i+1)-theta_list(2*i-1)-h/6*(omega_list(2*i-1)+4*omega_list(2*i)+omega_list(2*i+1));
        ceq8(i)=theta_list(2*i)-1/2*(theta_list(2*i-1)+theta_list(2*i+1))-h/8*(omega_list(2*i-1)-omega_list(2*i+1));
        ceq9(i)=omega_list(2*i+1)-omega_list(2*i-1)-h/6*(alpha_list(2*i-1)+4*alpha_list(2*i)+alpha_list(2*i+1));
        ceq10(i)=omega_list(2*i)-1/2*(omega_list(2*i-1)+omega_list(2*i+1))-h/8*(alpha_list(2*i-1)-alpha_list(2*i+1));
    end
else
    % Trapezoidal collocation
    h=T/(n-1); % time step size
    ceq3=x_list(2:end)-x_list(1:(end-1))-h/2*(v_list(1:(end-1))+v_list(2:end));
    ceq4=v_list(2:end)-v_list(1:(end-1))-h/2*(a_list(1:(end-1))+a_list(2:end));
    ceq5=theta_list(2:end)-theta_list(1:(end-1))-h/2*(omega_list(1:(end-1))+omega_list(2:end));
    ceq6=omega_list(2:end)-omega_list(1:(end-1))-h/2*(alpha_list(1:(end-1))+alpha_list(2:end));
    ceq7=0; ceq8=0; ceq9=0; ceq10=0;
end

ceq=[ceq1,ceq2,ceq3,ceq4,ceq5,ceq6,ceq7,ceq8,ceq9,ceq10];
% No inequal constraints c

end
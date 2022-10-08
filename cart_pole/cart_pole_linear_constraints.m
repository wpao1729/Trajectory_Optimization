function [Aeq,beq,A,b] = cart_pole_linear_constraints (n,n_dv,x1,xf,theta1,thetaf,T_fixed,T)

% Linear equality contraints ----------------------------------------------
n_bc=9; % number of boundary constraints
% n_cl=4*(n-1); % number of collocation constraints
Aeq=zeros(n_bc,n_dv); beq=zeros(n_bc,1); % initialize matrices

% Boundary constraints:
v1=0; % initial cart velocity (m/s)
vf=0; % final cart veloctiy (m/s)
omega1=0; % initial pole angular velocity (rad/s)
omegaf=0; % final pole angular velocity (rad/s)
Aeq(1,1)=1; Aeq(2,n)=1; Aeq(3,n+1)=1; Aeq(4,2*n)=1; Aeq(5,3*n+1)=1; Aeq(6,4*n)=1; Aeq(7,4*n+1)=1; Aeq(8,5*n)=1;
beq(1:8)=[x1,xf,v1,vf,theta1,thetaf,omega1,omegaf];
if T_fixed
    Aeq(9,7*n+1)=1; beq(9)=T;
end

% Linear inequality constraints -------------------------------------------
n_path=2*n; % number of path constraints
A=zeros(n_path,n_dv); b=zeros(n_path,1);

% Path constraints: force <= u_max
u_max = 100; % upper limit of the absolute value of force
A(1:n,(6*n+1):(7*n))=eye(n); A((n+1):(2*n),(6*n+1):(7*n))=-1*eye(n);
b(1:(2*n))=u_max;

end
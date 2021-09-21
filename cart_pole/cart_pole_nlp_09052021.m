clear all; clc;
tic;

mc=1; mp=1; L=0.5; g=9.81; %Parameters of the system

n=50; % number of points defining the path
t1=0; tf=2; t_list=linspace(t1,tf,n); h=(tf-t1)/(n-1); % Step size

% Decision Variables: x, v, a, theta, omega, alpha, u
n_bc=8;
Aeq=zeros(4*(n-1)+n_bc,7*n); beq=zeros(1,4*(n-1)+n_bc);

% Boundary constraints
x1=0; xf=1; v1=0; vf=0; theta1=-pi/2; thetaf=pi/2; omega1=0; omegaf=0;
Aeq(1,1)=1; Aeq(2,n)=1; Aeq(3,n+1)=1; Aeq(4,2*n)=1; Aeq(5,3*n+1)=1; Aeq(6,4*n)=1; Aeq(7,4*n+1)=1; Aeq(8,5*n)=1;
beq(1:n_bc)=[x1,xf,v1,vf,theta1,thetaf,omega1,omegaf];

% Collocation constraints - trapezoidal rule: f(k+1)-f(k)=h/2*(f'(k)+f'(k+1))
mat1=diag(-1*ones(1,n))+diag(ones(1,n-1),1); mat1=mat1(1:(end-1),:); % -1,1 matrix
mat2=diag(-h/2*ones(1,n))+diag(-h/2*ones(1,n-1),1); mat2=mat2(1:(end-1),:); %-h/2,-h/2 matrix
Aeq((n_bc+1):(n_bc+n-1),1:n)=mat1; Aeq((n_bc+1):(n_bc+n-1),(n+1):(2*n))=mat2; % x,v
Aeq((n_bc+n-1+1):(n_bc+2*(n-1)),(n+1):(2*n))=mat1; Aeq((n_bc+n):(n_bc+2*(n-1)),(2*n+1):(3*n))=mat2; % v,a
Aeq((n_bc+2*(n-1)+1):(n_bc+3*(n-1)),(3*n+1):(4*n))=mat1; Aeq((n_bc+2*(n-1)+1):(n_bc+3*(n-1)),(4*n+1):(5*n))=mat2; % theta,omega
Aeq((n_bc+3*(n-1)+1):(n_bc+4*(n-1)),(4*n+1):(5*n))=mat1; Aeq((n_bc+3*(n-1)+1):(n_bc+4*(n-1)),(5*n+1):(6*n))=mat2; % omega,alpha
    % beq for all collocation constraints equal zero.
    
% Path constraints: control
A=zeros(2*n,7*n);
A(1:n,(6*n+1):(7*n))=eye(n); A((n+1):(2*n),(6*n+1):(7*n))=-1*eye(n);
u_max=30; % upper limit of the absolute value of control
b=u_max*ones(1,2*n);

% Initialization
    a_x=-6*(xf-x1)/(tf^3);
    x_list_guess=a_x/3*(t_list.^3)-a_x*tf/2*(t_list.^2)+x1;
    v_list_guess=a_x*(t_list.^2)-a_x*tf*t_list;
    a_list_guess=2*a_x*t_list-a_x*tf;

    a_theta=-6*(thetaf-theta1)/(tf^3);
    theta_list_guess=a_theta/3*(t_list.^3)-a_theta*tf/2*(t_list.^2)+theta1;
    omega_list_guess=a_theta*(t_list.^2)-a_theta*tf*t_list;
    alpha_list_guess=2*a_theta*t_list-a_theta*tf;
    
    u_list_guess=(mc+mp)*a_list_guess-mp*L*sin(theta_list_guess).*alpha_list_guess-mp*L*cos(theta_list_guess).*(omega_list_guess.^2);
var_list_guess=[x_list_guess,v_list_guess,a_list_guess,theta_list_guess,omega_list_guess,alpha_list_guess,u_list_guess]; % Initial guess of the solution

% Optimization
options = optimoptions('fmincon','MaxFunctionEvaluations',200*n,'StepTolerance',1e-3);
func_cost=@(var_list)cart_pole_cost(var_list,n,tf);
func_nlcon=@(var_list)cart_pole_nlcon(var_list,n,mc,mp,L,g);
[var_list,cost]=fmincon(func_cost,var_list_guess,A,b,Aeq,beq,[],[],func_nlcon,options);
x_list=var_list(1:n); v_list=var_list((n+1):(2*n)); a_list=var_list((2*n+1):(3*n)); theta_list=var_list((3*n+1):(4*n)); omega_list=var_list((4*n+1):(5*n)); alpha_list=var_list((5*n+1):(6*n)); u_list=var_list((6*n+1):(7*n));

toc;

fprintf('\nT: %g.\n',tf);
fprintf('\nCost of the optimal path: %g.\n',cost);

figure(41);
L=0.5;
Cx_list=x_list; Cy_list=zeros(1,length(x_list)); Px_list=x_list+L*cos(theta_list); Py_list=L*sin(theta_list);
for k=1:length(x_list)
    plot([Cx_list(k),Px_list(k)],[Cy_list(k),Py_list(k)],'-bo');
    axis([-0.5,1.5,-1,1]);
    dt=tf/n; % Finish animation in tf
    pause(dt);
end


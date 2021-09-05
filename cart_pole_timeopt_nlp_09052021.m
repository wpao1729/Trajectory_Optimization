clear all; clc;
tic;

mc=1; mp=1; L=0.5; g=9.81; %Parameters of the system

n=200; % number of points defining the path

% Decision Variables: x, v, a, theta, omega, alpha, u, T
n_bc=8; n_dv=7*n+1;
Aeq=zeros(n_bc,n_dv); beq=zeros(1,n_bc);

% Boundary constraints
x1=0; xf=1; v1=0; vf=0; theta1=-pi/2; thetaf=pi/2; omega1=0; omegaf=0;
Aeq(1,1)=1; Aeq(2,n)=1; Aeq(3,n+1)=1; Aeq(4,2*n)=1; Aeq(5,3*n+1)=1; Aeq(6,4*n)=1; Aeq(7,4*n+1)=1; Aeq(8,5*n)=1;
beq(1:n_bc)=[x1,xf,v1,vf,theta1,thetaf,omega1,omegaf];
    
% Path constraints: control
A=zeros(2*n,n_dv);
A(1:n,(6*n+1):(7*n))=eye(n); A((n+1):(2*n),(6*n+1):(7*n))=-1*eye(n);
u_max=300; % upper limit of the absolute value of control
b=u_max*ones(1,2*n);

% Initialization
    T_guess=2; t_list_guess=linspace(0,T_guess,n);
    a_x=-6*(xf-x1)/(T_guess^3);
    x_list_guess=a_x/3*(t_list_guess.^3)-a_x*T_guess/2*(t_list_guess.^2)+x1;
    v_list_guess=a_x*(t_list_guess.^2)-a_x*T_guess*t_list_guess;
    a_list_guess=2*a_x*t_list_guess-a_x*T_guess;
    
    a_theta=-6*(thetaf-theta1)/(T_guess^3);
    theta_list_guess=a_theta/3*(t_list_guess.^3)-a_theta*T_guess/2*(t_list_guess.^2)+theta1;
    omega_list_guess=a_theta*(t_list_guess.^2)-a_theta*T_guess*t_list_guess;
    alpha_list_guess=2*a_theta*t_list_guess-a_theta*T_guess;
    
    u_list_guess=(mc+mp)*a_list_guess-mp*L*sin(theta_list_guess).*alpha_list_guess-mp*L*cos(theta_list_guess).*(omega_list_guess.^2);
var_list_guess=[x_list_guess,v_list_guess,a_list_guess,theta_list_guess,omega_list_guess,alpha_list_guess,u_list_guess,T_guess]; % Initial guess of the solution

% Optimization
options = optimoptions('fmincon','MaxFunctionEvaluations',200*n,'StepTolerance',1e-3);
func_cost=@(var_list)cart_pole_timeopt_cost(var_list,n);
func_nlcon=@(var_list)cart_pole_timeopt_nlcon(var_list,n,mc,mp,L,g);
[var_list,cost]=fmincon(func_cost,var_list_guess,A,b,Aeq,beq,[],[],func_nlcon,options);
x_list=var_list(1:n); v_list=var_list((n+1):(2*n)); a_list=var_list((2*n+1):(3*n)); 
theta_list=var_list((3*n+1):(4*n)); omega_list=var_list((4*n+1):(5*n)); alpha_list=var_list((5*n+1):(6*n)); 
u_list=var_list((6*n+1):(7*n)); T=var_list(7*n+1);
t1=0; tf=T; h=(tf-t1)/(n-1); t_list=linspace(t1,tf,n);

toc;

effort_total=h/2*(u_list(1)^2+u_list(end)^2+2*sum(u_list(2:(end-1)).^2));
fprintf('\nTotal effort: %g.\n',effort_total);
fprintf('\nT: %g.\n',T);
fprintf('\nCost of the optimal path (Weighted): %g.\n',cost);

figure(41);
Cx_list=x_list; Cy_list=zeros(1,length(x_list)); Px_list=x_list+L*cos(theta_list); Py_list=L*sin(theta_list);
for k=1:length(x_list)
    plot([Cx_list(k),Px_list(k)],[Cy_list(k),Py_list(k)],'-bo');
    axis([-0.5,1.5,-1,1]);
    dt=10/n; % Finish animation in 10s
    pause(dt);
end


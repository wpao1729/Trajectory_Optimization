clear all; clc;
tic;

mc=1; mp=1; L=0.5; g=9.81; %Parameters of the system

n=51; % number of points defining the path (must be odd)
n_seg=(n-1)/2; % number of segments. 1 segment contain 3 points: k, k+1/2, k+1
t1=0; tf=2; t_list=linspace(t1,tf,n); h=(tf-t1)/n_seg; % Step size (k to k+1)

% Decision Variables: x, v, a, theta, omega, alpha, u
n_bc=8;
Aeq=zeros((n_bc+8*n_seg),7*n); beq=zeros(1,(n_bc+8*n_seg));

% Boundary constraints
x1=0; xf=1; v1=0; vf=0; theta1=-pi/2; thetaf=pi/2; omega1=0; omegaf=0;
Aeq(1,1)=1; Aeq(2,n)=1; Aeq(3,n+1)=1; Aeq(4,2*n)=1; Aeq(5,3*n+1)=1; Aeq(6,4*n)=1; Aeq(7,4*n+1)=1; Aeq(8,5*n)=1;
beq(1:n_bc)=[x1,xf,v1,vf,theta1,thetaf,omega1,omegaf];

% Collocation constraints - Hermite-Simpson
mat1=zeros(n_seg,n); mat2=zeros(n_seg,n); mat3=zeros(n_seg,n); mat4=zeros(n_seg,n);
for i=1:n_seg
    mat1(i,(2*i-1))=-1; mat1(i,(2*i+1))=1;
    mat2(i,(2*i-1))=-h/6; mat2(i,(2*i))=-2/3*h; mat2(i,(2*i+1))=-h/6;
    mat3(i,(2*i-1))=-1/2; mat3(i,(2*i))=1; mat3(i,(2*i+1))=-1/2;
    mat4(i,(2*i-1))=-h/8; mat4(i,(2*i+1))=h/8;
end
mat0=[mat1,mat2;mat3,mat4];
Aeq((n_bc+1):(n_bc+2*n_seg),1:(2*n))=mat0; %x,v
Aeq((n_bc+2*n_seg+1):(n_bc+4*n_seg),(n+1):(3*n))=mat0; %v,a
Aeq((n_bc+4*n_seg+1):(n_bc+6*n_seg),(3*n+1):(5*n))=mat0; %theta,omega
Aeq((n_bc+6*n_seg+1):(n_bc+8*n_seg),(4*n+1):(6*n))=mat0; %omega,alpha
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
options = optimoptions('fmincon','MaxFunctionEvaluations',10000*n,'StepTolerance',1e-3);
func_cost=@(var_list)cart_pole_simpson_cost(var_list,n,h);
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


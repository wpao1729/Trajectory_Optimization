clear all; clc;

%% Parameters
n=30; % number of points defining the path
tf=1; % total time (s)
d=1; % distance the block moves (m)
m=1; % mass of the block (kg)

%% Formulate the trajectory optimization problem
% Create time list
t1=0; 
t_list=linspace(t1,tf,n); 
h=(tf-t1)/(n-1); % time step size

% Decision variables: x,v,a,u
n_dv=4*n; % number of decision variables

% Linear constraints
[Aeq,beq,A,b] = block_move_linear_constraints (n,n_dv,d,m,h);

% Cost function
func_cost=@(var_list)block_move_cost(var_list,tf,n);

% Initial guess
var_list_guess=[linspace(0,d,n),ones(1,n),ones(1,n),ones(1,n)]; % Initial guess of the solution

%% Solve optimization by fmincon
tic;
options = optimoptions('fmincon','Algorithm','interior-point','SubproblemAlgorithm','factorization',...
    'MaxFunctionEvaluations',2e6,'MaxIterations',1500,'StepTol',1e-8,'PlotFcn',{'optimplotfval','optimplotconstrviolation'}); % 'interior-point' 'factorization'
[var_list,cost]=fmincon(func_cost,var_list_guess,A,b,Aeq,beq,[],[],[],options);
toc;

%% Results Visualization
x_list_opt=var_list(1:n);
v_list_opt=var_list((n+1):2*n);
a_list_opt=var_list((2*n+1):3*n);
u_list_opt=var_list((3*n+1):4*n);

% True Analytic Solutions
tt = linspace(0,1,100);
x_list_ans=3*tt.^2-2*tt.^3;
v_list_ans=6*tt-6*tt.^2;
a_list_ans=6-12*tt;
u_list_ans=m*a_list_ans;

% Compare optimization results with analytic solutions
% Position x
figure(260);
plot(t_list,x_list_opt,'-bo',tt,x_list_ans,'-r');
grid minor;
xlabel('Time (s)'); ylabel('Position (m)');
legend('Optimal Path','True Analytic Solution');
title('Block Position x');

% Velocity v
figure(270);
plot(t_list,v_list_opt,'-bo',tt,v_list_ans,'-r');
grid minor;
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Optimal Path','True Analytic Solution');
title('Block Velocity v');

% Acceleration a
figure(280);
plot(t_list,a_list_opt,'-bo',tt,a_list_ans,'-r');
grid minor;
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
legend('Optimal Path','True Analytic Solution');
title('Block Acceleration a');

% Force u
figure(290);
plot(t_list,u_list_opt,'-bo',tt,u_list_ans,'-r');
grid minor;
xlabel('Time (s)'); ylabel('Force (N)');
legend('Optimal Path','True Analytic Solution');
title('Force u');
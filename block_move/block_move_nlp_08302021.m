clear all; clc;

n=30; % number of points defining the path
t1=0; tf=1; t_list=linspace(t1,tf,n); h=(tf-t1)/(n-1); % Step size

% Second approach: combine x_list and f_list, 14 deciding variables, 10 constraints
% xf_list=[x_list,f_list];
Aeq=zeros(n-1+4,2*n); beq=zeros(1,n-1+4);
% Boundary constraints
x1=0; xf=1; f1=0; ff=0;
Aeq(1,1)=1; Aeq(2,n)=1; Aeq(3,n+1)=1; Aeq(4,2*n)=1;
beq(1:4)=[x1,xf,f1,ff];
% Path constraints - trapezoidal rule: x(k+1)-x(k)=h/2*(f(k+1)+f(k))
mat1=diag(-1*ones(1,n))+diag(ones(1,n-1),1);
mat2=diag(-h/2*ones(1,n))+diag(-h/2*ones(1,n-1),1);
Aeq(5:end,1:n)=mat1(1:(end-1),:);
Aeq(5:end,(n+1):end)=mat2(1:(end-1),:);
beq(5:end)=zeros(1,n-1);
% Optimization by nonlinear programming fmincon
func=@(xf_list)compute_cost(t_list,xf_list);
xf_list_guess=[linspace(x1,xf,n),ones(1,n)]; % Initial guess of the solution
xf_list_opt=fmincon(func,xf_list_guess,[],[],Aeq,beq);

% Calculate the full representation of optimal path (state x, system dynamics f, control u)
x_list_opt=xf_list_opt(1:n);
f_list_opt=xf_list_opt((n+1):end);
u_list_opt=(f_list_opt(2:end)-f_list_opt(1:(end-1)))/h;
cost_opt=h*sum(u_list_opt.^2);

% True Analytic Solutions
tt = linspace(0,1,100);
x_list_ans=3*tt.^2-2*tt.^3;
f_list_ans=6*tt-6*tt.^2;
u_list_ans=6-12*tt;

% Visualize optimization results
figure(26);
spline1=spline(t_list,[f1,x_list_opt,ff]);
plot(t_list,x_list_opt,'bo',tt,ppval(spline1,tt),'-b',tt,x_list_ans,'-r');
grid minor;
xlabel('Time'); ylabel('Position');
legend('Optimal Path','Spline Interpolation','True Analytic Solution');
title('State x');

figure(27);
plot(t_list,f_list_opt,'-bo',tt,f_list_ans,'-r');
grid minor;
xlabel('Time'); ylabel('Velocity');
legend('Optimal Path','True Analytic Solution');
title('System Dynamics f');

figure(28);
t_list_2=(t_list(1:(end-1))+t_list(2:end))/2;
plot(t_list_2,u_list_opt,'-bo',tt,u_list_ans,'-r');
grid minor;
xlabel('Time'); ylabel('Force');
legend('Optimal Path','True Analytic Solution');
title('Control u');
clear all; clc;

n=5; % number of intermediate points
t1=0; tf=1; t_list=linspace(t1,tf,n+2); h=(tf-t1)/(n+1);
x1=0; xf=1;
fi=0; ff=0;

x0=1/6.*(1:(n-1)); % Initial guess of the solution
% Optimization by nonlinear programming fmincon---------------------------
x_inter_opt=fmincon(@compute_cost,x0,[],[],[],[],zeros(1,4),ones(1,4));
% ------------------------------------------------------------------------

% Calculate the full representation of optimal path (state x, system dynamics f, control u)
x_list_opt=[x1,x_inter_opt,-1,xf]; % x(end-1) unassigned
f_list_opt=[fi,zeros(1,n),ff];
for k=1:(n-1)
   f_list_opt(k+1)=2/h*(x_list_opt(k+1)-x_list_opt(k))-f_list_opt(k);
end
f_list_opt(end-1)=(x_list_opt(end)-x_list_opt(end-2))/h-1/2*(f_list_opt(end)+f_list_opt(end-2));
x_list_opt(end-1)=x_list_opt(end)-1/2*h*(f_list_opt(end-1)+f_list_opt(end));
u_list_opt=(f_list_opt(2:end)-f_list_opt(1:(end-1)))/h;

tt = linspace(0,1,100);
x_list_ans=3*tt.^2-2*tt.^3;
f_list_ans=6*tt-6*tt.^2;
u_list_ans=6-12*tt;

figure(26);
spline1=spline(t_list,[0 x_list_opt 0]);
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

function cost = compute_cost(x_inter)
    x_list=[x1,x_inter,-1,xf]; % x(end-1) unassigned
    f_list=[fi,zeros(1,n),ff];
    for k=1:(n-1)
       f_list(k+1)=2/h*(x_list(k+1)-x_list(k))-f_list(k);
    end
    f_list(end-1)=(x_list(end)-x_list(end-2))/h-1/2*(f_list(end)+f_list(end-2));
    x_list(end-1)=x_list(end)-1/2*h*(f_list(end-1)+f_list(end));
    u_list=(f_list(2:end)-f_list(1:(end-1)))/h;
    cost=h*sum(u_list.^2);
end
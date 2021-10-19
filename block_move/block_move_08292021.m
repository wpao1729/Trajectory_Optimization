% Trajectory Optimization. Block-move example.
% System Dynamics: x'=f
% Trapezoidal Rule: x(k+1)-x(k)=1/2*h*(f(k)+f(k+1))
% f(k+1)-f(k)=1/2*h*(u(k)+u(k+1))
% Boundary constraints: x(0)=0, f(0)=0; x(1)=1, f(1)=0.
% Path constraints: x(k+1)>=x(k), f(:)>=0;
% Minimum-force objective: min(integral(u^2))

n=5; % number of intermediate discrete points
x1=0; xf=1; 
h=(xf-x1)/(n+1); % time step
% Boundary constraints
t_list=linspace(0,1,n+2);
x_list(1)=x1; x_list(n+2)=xf;
f_list(1)=0; f_list(n+2)=0;

accuracy=0.005;
cost_min=1000;
for x2=x1:accuracy:xf
    for x3=x2:accuracy:xf
        for x4=x3:accuracy:xf
            for x5=x4:accuracy:xf
                % Trial --------------------------------------------------
                % Construct vector of state x,f
                x6=0;
                x_list=[x1,x2,x3,x4,x5,x6,xf];
                f_list(1)=0; f_list(n+2)=0;
                for k=1:(length(x_list)-3)
                   f_list(k+1)=2/h*(x_list(k+1)-x_list(k))-f_list(k);
                end
                f_list(end-1)=(x_list(end)-x_list(end-2))/h-1/2*(f_list(end)+f_list(end-2));
                x_list(end-1)=x_list(end)-1/2*h*(f_list(end-1)+f_list(end));
                % Check if all f>=0 to satisfy path constraint
                if any(f_list<0) 
                    continue
                end
                % Construct vector of control u
%                 u1=((f_list(2)-f_list(1))/h);
%                 u_list=zeros(1,length(x_list)); u_list(1)=u1;
%                 for i=1:(length(u_list)-1)
%                     u_list(i+1)=2/h*(f_list(i+1)-f_list(i))-u_list(i);
%                 end
                u_list=(f_list(2:end)-f_list(1:(end-1)))/h; % length(u_list)=length(x_list)-1
 
                % Find minimum-force objective
                % cost=1/2*h*(u_list(1)^2+2*sum(u_list(2:end-1).^2)+u_list(end)^2); % cost=integral(u^2)
                % cost=1/3*h*(u_list(1)^2+2*sum(u_list(2:end-1).^2)+u_list(end)^2+sum(u_list(1:end-1).*u_list(2:end)));
                cost=h*sum(u_list.^2);
                if cost<cost_min
                    x_list_opt=x_list; f_list_opt=f_list; u_list_opt=u_list;
                    cost_min=cost;
                end
                % --------------------------------------------------------
            end
        end
    end
end

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
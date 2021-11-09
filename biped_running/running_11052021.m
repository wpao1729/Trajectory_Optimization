clear all; clc;

%% Parameters
% User defined parameters:  
tf=0.6; % time for one step (s)
d_run=0.8; % distance of one running step from the previous stance (m)
h_run=0; % elevation of one running step from the previous stance (m)
miu=0.8; % friction coefficient of the terrain
ft_clr=0.005; % foot clearance (m)
iswalking=true; % set true if walking, false if running
n=15; % Number of knot points including 2 ends

% Physics parameters:
g=9.81; % gravitational acceleration (m/s^2)
m_list=[3.2,6.8,20,6.8,3.2]; % mass of links (kg)
I_list=[0.93,1.08,2.22,1.08,0.93]; % moment of inertia of links (kg*m^2)
L_list=[0.4,0.4,0.625,0.4,0.4]; % Length of links (m)
d_list=[0.128,0.163,0.2,0.163,0.128]; % distance of link (center of mass) to the nearest joint (m)

%% Formulate the trajectory optimization problem
% Decision variables:
% q(5*n),vq(5*n),aq(5*n),P(12*n),vP0(2*n),aP0(2*n),vP5(2*n),aP5(2*n),G(10*n),vG(10*n),aG(10*n),GRF(4*n),u(6*n),t1,t2
n_dv=75*n+2; % number of decision variables

% Create time array
t1=0;
t_list=linspace(t1,tf,n); % time array
h=(tf-t1)/(n-1); % time step size

% Linear Constraints:
[Aeq,beq,A,b] = running_linear_constraints(n,m_list,g,d_run,h_run,tf,h,n_dv,iswalking,miu);

% Nonlinear Constraints:
func_nlcon=@(var_list)running_nonlinear_constraints(var_list,n,m_list,I_list,g,L_list,d_list,d_run,h_run,tf,h,ft_clr);

% Cost function:
func_cost=@(var_list)running_cost(var_list,n,tf);

% Initialization
qi_list=[20,-20,0,-15,45]/180*pi; % initial guess for initial joint angles
qf_list=fliplr(qi_list); % intial guess for  final joint angles
var_list_guess=initialize_running(n,qi_list,qf_list,L_list,d_list,m_list,I_list,g,d_run,h_run,t_list,iswalking);

%% Solve the trajectory optimization problem with fmincon
tic;
options = optimoptions('fmincon','Algorithm','interior-point','SubproblemAlgorithm','cg',...
    'MaxFunctionEvaluations',1e6,'MaxIterations',1000,'StepTol',1e-8,'TolFun',1e-3,'TolCon',1e-4,'PlotFcn',{'optimplotfval','optimplotconstrviolation'}); 
[var_list,cost]=fmincon(func_cost,var_list_guess,A,b,Aeq,beq,[],[],func_nlcon,options);
toc;

%% Animation
n_step=5; % Number of steps for animation

% Create robot model shape for animation
P_list=var_list((15*n+1):(27*n)); 
P_mat=reshape(P_list,12,n);
x_mat=zeros(7,n_step*(n-1));
y_mat=zeros(7,n_step*(n-1));
for i=1:n_step
    x_mat(:,((i-1)*(n-1)+1):(i*(n-1)))=[P_mat([1,3,5,7,5,9,11],1:(n-1))]+(i-1)*d_run;
    y_mat(:,((i-1)*(n-1)+1):(i*(n-1)))=[P_mat([2,4,6,8,6,10,12],1:(n-1))]+(i-1)*h_run;
end

% Animate
fig_animation=figure(483); clf(fig_animation);
clear frame_list;
but=uicontrol(fig_animation,'Style','pushbutton','Position',[10 10 50 20],'String','Play','Callback',@(source,event)play(source,event,fig_animation,x_mat,y_mat,n_step,d_run,h_run));
for k=1:size(x_mat,2)
    plot(x_mat(:,k),y_mat(:,k),'-bo'); hold on;
    for kk=1:n_step
    plot([kk*d_run-0.1,kk*d_run-0.1,(kk+1)*d_run-0.1],[(kk-1)*h_run,kk*h_run,kk*h_run]); 
    end
    hold off;
    axis([-0.5,(n_step+1)*d_run,0,1.5+(n_step)*h_run]); set(gca,'DataAspectRatio',[1 1 1]);
    title('Running Animation','FontSize', 15);
    xlabel('x (m)','FontSize', 15); ylabel('z (m)','FontSize', 15);
    frame_list(k) = getframe(gcf) ;
    drawnow;
end

%% Save animation to MP4
folder='C:/Users/5baow/OneDrive - Georgia Institute of Technology/Desktop/LIDAR Gatech/MATLAB Workspace/Matlab Results/Running Animation/';
if iswalking
    filename=sprintf('walking_d%gm_h%gm_t%gs_cost%g.mp4',d_run,h_run,tf,cost);
else
    filename=sprintf('running_d%gm_h%gm_t%gs_cost%g.mp4',d_run,h_run,tf,cost);
end
writerObj = VideoWriter([folder,filename],'MPEG-4'); % for mp4 file
framerate=1/(tf/n); % frame rate (frequency)
writerObj.FrameRate = framerate;
open(writerObj);
for i=1:length(frame_list) 
    writeVideo(writerObj, frame_list(i));
end
close(writerObj);

%% Plot torques vs time
u_mat=reshape(var_list((69*n+1):(75*n)),6,n);
fig_torque=figure(512);
plot(t_list,u_mat(1,:)); hold on; plot(t_list,u_mat(2,:)); plot(t_list,u_mat(3,:)); plot(t_list,u_mat(4,:)); plot(t_list,u_mat(5,:)); plot(t_list,u_mat(6,:)); hold off;
xlabel('time(s)'); ylabel('Torque (N*m)');
legend('u1','u2','u3','u4','u5','u6');
title('Torques');

%% Functions section
function play(source,event,fig,x_mat,y_mat,n_step,d_run,h_run)
    fig;
    for k=1:size(x_mat,2)
        plot(x_mat(:,k),y_mat(:,k),'-bo'); hold on;
        for kk=1:n_step
        plot([kk*d_run-0.1,kk*d_run-0.1,(kk+1)*d_run-0.1],[(kk-1)*h_run,kk*h_run,kk*h_run]); 
        end
        hold off;
        axis([-0.5,(n_step+1)*d_run,0,1.5+(n_step)*h_run]); set(gca,'DataAspectRatio',[1 1 1]);
        title('Running Animation','FontSize', 15);
        xlabel('x (m)','FontSize', 15); ylabel('z (m)','FontSize', 15);
        drawnow;
    end
end

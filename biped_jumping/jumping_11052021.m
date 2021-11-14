clear; clc;

%% Parameters
% User-defined parameters
tf=2; % total time (s) % tf~=2s for flip
d_jump=0.5; % distance of jump (m)
h_jump=0; % elevation of jump (m)
miu=0.8; % friction coefficient of the terrain
n=20; % Number of knot points including 2 ends
ft_clr=0.005; % foot clearance (m)
qi_list=[0,0,0]/180*pi; % starting position
isflip=true; % set true if flip, false if normal jump

% Physics parameters:
g=9.81; % gravitational acceleration (m/s^2)
m_list=[3.2*2,6.8*2,20]; % mass of links (kg)
I_list=[0.93*2,1.08*2,2.22]; % moment of inertia of links (kg*m^2)
L_list=[0.4,0.4,0.625]; % length of links (m)
d_list=[0.272,0.237,0.2]; % distance of link (center of mass) to the lower joint (m)

%% Formulate the trajectory optimization problem
% Decision variables:
% q(3*n),vq(3*n),aq(3*n),P(4*2*n),vP0(2*n),aP0(2*n),G(3*2*n),vG(3*2*n),aG(3*2*n),u(3*n),t1,t2
n_dv=42*n+2; % number of decision variables

% Create time array
t1=0;
t_list=linspace(t1,tf,n); % time array
h=(tf-t1)/(n-1); % time step size

% Create final position of q123
if isflip
    qf_list=[360,360,360]/180*pi; % final q123 for front flip
else
    qf_list=[0,0,0]/180*pi; % final q123 for normal jumping
end

% Linear constraints
[Aeq,beq,A,b] = jumping_linear_constraints (n,d_jump,h_jump,tf,h,qi_list,qf_list,n_dv);

% Nonlinear constraints
func_nlcon=@(var_list)jumping_nonlinear_constraints(var_list,n,m_list,I_list,g,L_list,d_list,d_jump,h_jump,miu,tf,h,ft_clr,isflip);

% Cost function
func_cost=@(var_list)jumping_cost(var_list,n,tf);

% Initialization
var_list_guess=initialize_jumping(n,qi_list,qf_list,L_list,d_list,m_list,I_list,g,d_jump,h_jump,t_list,isflip);

%% Solve optimization with fmincon
tic;
options = optimoptions('fmincon','Algorithm','interior-point','SubproblemAlgorithm','cg',...
    'MaxFunctionEvaluations',2e6,'MaxIterations',1500,'StepTol',1e-8,'PlotFcn',{'optimplotfval','optimplotconstrviolation'}); % 'interior-point' 'factorization'
[var_list,cost]=fmincon(func_cost,var_list_guess,A,b,Aeq,beq,[],[],func_nlcon,options);
toc;

%% Animation
P_list=var_list((9*n+1):(17*n)); 
P_mat=reshape(P_list,8,n);

fig_animation=figure(614); clf(fig_animation);
clear frame_list;
but=uicontrol(fig_animation,'Style','pushbutton','Position',[10 10 50 20],'String','Play','Callback',@(source,event)play(source,event,fig_animation,P_mat,d_jump,h_jump,n));
for k=1:n
    plot(P_mat(1:2:end,k)',P_mat(2:2:end,k)','-bo'); hold on;
    plot([-0.5,d_jump-0.1,d_jump-0.1,d_jump+0.5],[0,0,h_jump,h_jump],'-r'); hold off;
    %plot([d_jump-0.5,d_jump+0.1,d_jump+0.1,0],[h_jump,h_jump,0,0]); hold off;
    axis([-1,d_jump+1,0,h_jump+2]); 
    set(gca,'DataAspectRatio',[1 1 1]);
    title('Jumping Animation','FontSize', 15);
    xlabel('x (m)','FontSize', 15); ylabel('z (m)','FontSize', 15);
    frame_list(k) = getframe(gcf) ;
    drawnow;
end

%% Save animation to mp4 file
% Specify folder to save the animation
browse_for_folder=true; % set true to browse for folder, false to manually specify
if browse_for_folder
    export_folder=uigetdir; % open folder selection box
else
    % User manually specify the destination folder
    workspace=fileparts(fileparts(pwd));
    export_folder=fullfile(workspace,'Matlab Results','Jumping Animation');
end
fprintf('\nAnimation Export folder: %s\n',export_folder);

% Specify animation file name
if isflip
    filename=sprintf('frontflip_d%gm_h%gm_t%gs_cost%g.mp4',d_jump,h_jump,tf,cost);
else
    filename=sprintf('jumping_d%gm_h%gm_t%gs_cost%g.mp4',d_jump,h_jump,tf,cost);
end

writerObj = VideoWriter(fullfile(export_folder,filename),'MPEG-4'); % for mp4 file
framerate=1/(tf/n); % frame rate (frequency)
writerObj.FrameRate = framerate;
open(writerObj);
for i=1:length(frame_list) 
    writeVideo(writerObj, frame_list(i));
end
close(writerObj);

%% Plot torques vs time
u_mat=reshape(var_list((39*n+1):(42*n)),3,n);
figure(63);
plot(t_list,u_mat(1,:)); hold on; plot(t_list,u_mat(2,:)); plot(t_list,u_mat(3,:)); hold off
grid minor;
xlabel('time(s)'); ylabel('Torque');
legend('u1','u2','u3');

%% Functions section
function play(source,event,fig,P_mat,d_jump,h_jump,n)
    fig;
    for k=1:n
        plot(P_mat(1:2:end,k)',P_mat(2:2:end,k)','-bo'); hold on;
        plot([-0.5,d_jump-0.1,d_jump-0.1,d_jump+0.5],[0,0,h_jump,h_jump],'-r'); hold off;
        %plot([d_jump-0.5,d_jump+0.1,d_jump+0.1,0],[h_jump,h_jump,0,0]); hold off;
        axis([-1,d_jump+1,0,h_jump+2]); 
        set(gca,'DataAspectRatio',[1 1 1]);
        title('Jumping Animation','FontSize', 15);
        xlabel('x (m)','FontSize', 15); ylabel('z (m)','FontSize', 15);
        drawnow;
    end
end

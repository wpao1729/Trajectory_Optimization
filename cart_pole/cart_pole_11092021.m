clear; clc;

%% Parameters
% User-defined parameters:
n=50; % number of knot (include 2 ends)
d=1; % distance the cart moves (m)
T=3; % total time (s)
T_fixed=true; % set true if T is fixed, false if optimize T
simpson=true; % set true if collocation by simpson, false if collocation by trapezoidal

% Physics parameters:
m_cart=1; % mass of cart (kg)
m_pole=1; % mass of pole (kg)
L=0.5; % length of the rod connecting the cart and the pole (m)
g=9.81; % gravitational acceleration (m/s^2)

%% Formulate the trajectory optimization problem
% Simpson collocation check: n has to be odd number.
if simpson
    if mod(n,2)==0
        n=n+1;
        fprintf('\nNumber of knot points is updated to %g for simpson collocation.\n',n);
    end
end

% Decision Variables: x, v, a, theta, omega, alpha, u, T
n_dv=7*n+1; % number of decision variables

% Cart-pole problem:
x1=0; % initial cart position (m)
xf=d; % final cart position (m)
theta1=-pi/2; % initial pole angle (rad)
thetaf=pi/2; % final pole angle (rad)

% Linear constraints:
[Aeq,beq,A,b] = cart_pole_linear_constraints (n,n_dv,x1,xf,theta1,thetaf,T_fixed,T);

% Nonlinear constraints:
func_nlcon=@(var_list)cart_pole_nonlinear_constraints(var_list,n,m_cart,m_pole,L,g,simpson);

% Cost function:
func_cost=@(var_list)cart_pole_cost(var_list,n,simpson);

% Initialization
var_list_guess = initialize_cart_pole (n,x1,xf,theta1,thetaf,m_cart,m_pole,L);

%% Solve the optimization problem by fmincon
tic;
options = optimoptions('fmincon','Algorithm','interior-point','SubproblemAlgorithm','factorization',...
    'MaxFunctionEvaluations',2e6,'MaxIterations',1500,'StepTol',1e-8,'PlotFcn',{'optimplotfval','optimplotconstrviolation'}); % 'interior-point' 'factorization'
[var_list,cost]=fmincon(func_cost,var_list_guess,A,b,Aeq,beq,[],[],func_nlcon,options);
toc;

%% Animation
x_list=var_list(1:n);
theta_list=var_list((3*n+1):(4*n));
u_list=var_list((6*n+1):(7*n)); 
T=var_list(7*n+1);
h=(T-0)/(n-1); t_list=linspace(0,T,n);

% Animation
Cx_list=x_list; Cy_list=zeros(1,length(x_list)); Px_list=x_list+L*cos(theta_list); Py_list=L*sin(theta_list);
l_cart=0.15; h_cart=0.1; 
fig_animation=figure(411); clf(fig_animation);
clear frame_list;
but=uicontrol(fig_animation,'Style','pushbutton','Position',[10 10 50 20],'String','Replay','Interruptible','off','BusyAction','cancel','Callback',@(source,event)play(source,event,fig_animation,Cx_list,Cy_list,Px_list,Py_list,l_cart,h_cart,n,d,L));
for k=1:n
    plot([Cx_list(k),Px_list(k)],[Cy_list(k),Py_list(k)],'-bo'); hold on;
    plot([Cx_list(k)-l_cart/2,Cx_list(k)+l_cart/2,Cx_list(k)+l_cart/2,Cx_list(k)-l_cart/2,Cx_list(k)-l_cart/2],...
        [Cy_list(k)-h_cart/2,Cy_list(k)-h_cart/2,Cy_list(k)+h_cart/2,Cy_list(k)+h_cart/2,Cy_list(k)-h_cart/2],'-r');
    plot([-1,d+1],[0,0],'-k'); hold off;
    axis([-1,d+0.5,-L-0.2,L+0.2]);
    set(gca,'DataAspectRatio',[1 1 1]);
    title('Cart-pole Animation','FontSize', 15);
    xlabel('x (m)','FontSize', 15); ylabel('z (m)','FontSize', 15);
    frame_list(k) = getframe(gcf) ;
    drawnow;
end

%% Save animation to mp4
% Specify folder to save the animation
browse_for_folder=true; % set true to browse for folder, false to manually specify
if browse_for_folder
    export_folder=uigetdir; % open folder selection box
else
    % User manually specify the destination folder
    workspace=fileparts(fileparts(pwd));
    export_folder=fullfile(workspace,'Matlab Results','Cart-pole Animation');
end
fprintf('\nAnimation Export folder: %s\n',export_folder);

filename=sprintf('cart_pole_d%gm_t%gs_cost%g.mp4',d,T,cost);

writerObj = VideoWriter(fullfile(export_folder,filename),'MPEG-4'); % for mp4 file
% writerObj = VideoWriter([folder,filename],'MPEG-4'); % for mp4 file
framerate=1/(T/n); % frame rate (frequency)
writerObj.FrameRate = framerate;
open(writerObj);
for i=1:length(frame_list) 
    writeVideo(writerObj, frame_list(i));
end
close(writerObj);

%% Plot force u over time
u_list=var_list((6*n+1):(7*n)); 
figure(209);
plot(t_list,u_list);
grid minor;
xlabel('time(s)'); ylabel('Force(N)');
legend('u');

% % Print
% effort_total=h/2*(u_list(1)^2+u_list(end)^2+2*sum(u_list(2:(end-1)).^2));
% fprintf('\nTotal effort: %g.\n',effort_total);
% fprintf('\nT: %g.\n',T);
% fprintf('\nCost of the optimal path (Weighted): %g.\n',cost);

%% Functions section
function play(source,event,fig,Cx_list,Cy_list,Px_list,Py_list,l_cart,h_cart,n,d,L)
    fig;
    for k=1:n
        plot([Cx_list(k),Px_list(k)],[Cy_list(k),Py_list(k)],'-bo'); hold on;
        plot([Cx_list(k)-l_cart/2,Cx_list(k)+l_cart/2,Cx_list(k)+l_cart/2,Cx_list(k)-l_cart/2,Cx_list(k)-l_cart/2],...
            [Cy_list(k)-h_cart/2,Cy_list(k)-h_cart/2,Cy_list(k)+h_cart/2,Cy_list(k)+h_cart/2,Cy_list(k)-h_cart/2],'-r');
        plot([-1,d+1],[0,0],'-k'); hold off;
        axis([-1,d+0.5,-L-0.2,L+0.2]);
        set(gca,'DataAspectRatio',[1 1 1]);
        title('Cart-pole Animation','FontSize', 15);
        xlabel('x (m)','FontSize', 15); ylabel('z (m)','FontSize', 15);
        drawnow;
    end
end

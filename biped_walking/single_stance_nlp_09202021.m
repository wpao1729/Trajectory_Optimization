clear all; clc;

tic;

% Decision variables:
% q(1,2,3,4,5),vq(1,2,3,4,5),aq(1,2,3,4,5),P(i,j)(1,2,3,4,5),G(i,j)(1,2,3,4,5),vG(i,j)(1,2,3,4,5),aG(i,j)(1,2,3,4,5),u(1,2,3,4,5)

n=6; % Number of knot points including 2 ends
t1=0; tf=1; t_list=linspace(t1,tf,n); h=(tf-t1)/(n-1); % Step size
qi_list=[-10,-20,5,10,20]/180*pi; qf_list=[20,10,5,-20,-10]/180*pi;

% Model parameters:
g=9.81;
% m_list=[1,3,10,3,1]; L_list=[0.5,0.5,0.5,0.5,0.5]; d_list=[0.25,0.25,0.25,0.25,0.25]; I_list=1/12*m_list.*(L_list.^2);
m_list=[3.2,6.8,20,6.8,3.2]; I_list=[0.93,1.08,2.22,1.08,0.93]; L_list=[0.4,0.4,0.625,0.4,0.4]; d_list=[0.272,0.237,0.2,0.163,0.128];

% Linear constraints------------------------------------------------------
n_bc=20;
Aeq=zeros(30*(n-1)+n_bc+n,60*n); beq=zeros(30*(n-1)+n_bc+n,1);

% Collocation constraints
mat1=diag(-1*ones(1,5*n))+diag(ones(1,5*(n-1)),5); mat1=mat1(1:(end-5),:); % -1,1 matrix, 5*(n-1) x 5*n
mat2=diag(-h/2*ones(1,5*n))+diag(-h/2*ones(1,5*(n-1)),5); mat2=mat2(1:(end-5),:); %-h/2,-h/2 matrix, 5*(n-1) x 5*n
mat11=diag(-1*ones(1,10*n))+diag(ones(1,10*(n-1)),10); mat11=mat11(1:(end-10),:); % -1,1 matrix, 10*(n-1) x 10*n
mat22=diag(-h/2*ones(1,10*n))+diag(-h/2*ones(1,10*(n-1)),10); mat22=mat22(1:(end-10),:); %-h/2,-h/2 matrix, 10*(n-1) x 10*n
Aeq(1:(5*(n-1)),1:(5*n))=mat1; Aeq(1:(5*(n-1)),(5*n+1):(10*n))=mat2; % q12345 % vq12345
Aeq((5*(n-1)+1):(10*(n-1)),(5*n+1):(10*n))=mat1; Aeq((5*(n-1)+1):(10*(n-1)),(10*n+1):(15*n))=mat2; % vq12345 & aq12345
% Aeq((10*(n-1)+1):(20*(n-1)),(25*n+1):(35*n))=mat11; Aeq((10*(n-1)+1):(20*(n-1)),(35*n+1):(45*n))=mat22; % Gij12345 & vGij12345
% Aeq((20*(n-1)+1):(30*(n-1)),(35*n+1):(45*n))=mat11; Aeq((20*(n-1)+1):(30*(n-1)),(45*n+1):(55*n))=mat22; % vGij12435 & aGij12345
    %beq for collocation constraints are all 0.

% Boundary constraints
Aeq_bc=zeros(n_bc,60*n); beq_bc=zeros(n_bc,1);
Aeq_bc(1:5,1:5)=eye(5); beq_bc(1:5)=qi_list'; 
Aeq_bc(6:10,(5*n-4):(5*n))=eye(5); beq_bc(6:10)=qf_list';
Aeq_bc(11:15,(5*n+1):(5*n+5))=eye(5); %beq=0;
Aeq_bc(16:20,(10*n-4):(10*n))=eye(5); %beq=0;
Aeq((30*(n-1)+1):30*(n-1)+n_bc,:)=Aeq_bc; beq((30*(n-1)+1):30*(n-1)+n_bc)=beq_bc;

% Special constraints: u1=0;
Aeq_bc((30*(n-1)+n_bc+1):(30*(n-1)+n_bc+n),(55*n+1):5:(60*n))=eye(n);
    %beq for special constraints are all 0;

% Inequal path constraints: P5j>=0, -P5j<=0
h_clearance=0.02;
A=zeros(n-2,60*n); b=zeros(n-2,1);
A(:,(15*n+20):10:(25*n-10))=-1*eye(n-2); b(:)=-h_clearance;

% ------------------------------------------------------------------------

% Initialization
var_list_guess=initialize(n,t_list,qi_list,qf_list,L_list,d_list,m_list,I_list,g);

% !!! Use get initialization with linspace, dynamics and kinematics !!!

% Optimization with fmincon
options = optimoptions('fmincon','MaxFunctionEvaluations',10e4,'StepTol',1e-6);
func_cost=@(var_list)single_stance_cost(var_list,n,tf);
func_nlcon=@(var_list)single_stance_nlcon(var_list,n,m_list,I_list,g,L_list,d_list);
[var_list,cost]=fmincon(func_cost,var_list_guess,A,b,Aeq,beq,[],[],func_nlcon,options);

toc;

fprintf('\nCost of the optimal path: %g.\n',cost);

% % Visualize
P_list=var_list((15*n+1):(25*n)); 
P_mat=reshape(P_list,10,n);
% figure(34);
% for k=1:n
%     plot([0,P_mat([1,3,5,3,7,9],k)'],[0,P_mat([2,4,6,4,8,10],k)'],'-bo');
%     axis([-1,1,0,1.6]);
%     dt=tf/n; % Finish animation in tf
%     pause(dt);
% end

figure(35);
for k=1:n
    plot([0,P_mat([1,3,5,3,7,9],k)'],[0,P_mat([2,4,6,4,8,10],k)'],'-o'); hold on;
    axis([-0.5,0.5,0,1.5]);
    dt=tf/n; % Finish animation in tf
    pause(dt);
end
hold off;
legend('t=0s','t=0.2s','t=0.4s','t=0.6s','t=0.8s','t=1s','FontSize', 14);
title('Single-stance Walking in 1 sec','FontSize', 18);
xlabel('x (m)','FontSize', 18); ylabel('z (m)','FontSize', 18);

u_mat=reshape(var_list((55*n+1):(60*n)),5,n);
effort_list=sum(u_mat.^2);
figure(36);
plot(t_list,effort_list,'-r');
xlabel('time(s)'); ylabel('Sum torque squared');




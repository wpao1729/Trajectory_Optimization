clear all; clc;

tic;

t1=0; tf=1; d_jump=1; h_jump=0.2;

n=20; % Number of knot points including 2 ends
% Decision variables:
% q(3*n),vq(3*n),aq(3*n),P(4*2*n),vP0(2*n),aP0(2*n),G(3*2*n),vG(3*2*n),aG(3*2*n),u(3*n),t1,t2
n_dv=42*n; % (+2 time variables for phase transition)

t_list=linspace(t1,tf,n); h=(tf-t1)/(n-1); % Step size
qi_list=[0,0,0]/180*pi; qf_list=qi_list;
miu=0.8; % coefficient of friction

% Model parameters:
g=9.81;
m_list=[3.2*2,6.8*2,20]; I_list=[0.93*2,1.08*2,2.22]; L_list=[0.4,0.4,0.625]; d_list=[0.272,0.237,0.2];
m1=m_list(1); m2=m_list(2); m3=m_list(3);

% Linear constraints------------------------------------------------------
n_cl=6*(n-1); % collocation between (q,vq) and (vq,aq)
n_bc=30; % 2 for P5ij(n)=(D,0), 2 for P5ij(n-1)=P5ij(n)
% n_rp=0; % repeatability
Aeq=zeros(n_cl+n_bc,n_dv); beq=zeros(n_cl+n_bc,1);

% Collocation constraints (trapezoidal)
mat1=diag(-1*ones(1,3*n))+diag(ones(1,3*(n-1)),3); mat1=mat1(1:(end-3),:); % -1,1 matrix, 3*(n-1) x 3*n
mat2=diag(-h/2*ones(1,3*n))+diag(-h/2*ones(1,3*(n-1)),3); mat2=mat2(1:(end-3),:); %-h/2,-h/2 matrix, 3*(n-1) x 3*n
Aeq(1:(3*(n-1)),1:(3*n))=mat1; Aeq(1:(3*(n-1)),(3*n+1):(6*n))=mat2; % q123 % vq123
Aeq((3*(n-1)+1):(6*(n-1)),(3*n+1):(6*n))=mat1; Aeq((3*(n-1)+1):(6*(n-1)),(6*n+1):(9*n))=mat2; % vq123 & aq123
% Collocation on center of mass
% mat3=zeros(2*n,6*n);
% mat3(1:2:end,1:6:end)=diag(-m1*ones(1,n))+diag(m1*ones(1,n-1),1); mat3(1:2:end,3:6:end)=diag(-m2*ones(1,n))+diag(m2*ones(1,n-1),1); mat3(1:2:end,5:6:end)=diag(-m3*ones(1,n))+diag(m3*ones(1,n-1),1);
% mat3(2:2:end,2:6:end)=diag(-m1*ones(1,n))+diag(m1*ones(1,n-1),1); mat3(2:2:end,4:6:end)=diag(-m2*ones(1,n))+diag(m2*ones(1,n-1),1); mat3(2:2:end,6:6:end)=diag(-m3*ones(1,n))+diag(m3*ones(1,n-1),1);
% mat3=mat3(1:(end-2),:);
% mat4=zeros(2*n,6*n);
% mat4(1:2:end,1:6:end)=diag(-h/2*m1*ones(1,n))+diag(-h/2*m1*ones(1,n-1),1); mat4(1:2:end,3:6:end)=diag(-h/2*m2*ones(1,n))+diag(-h/2*m2*ones(1,n-1),1); mat4(1:2:end,5:6:end)=diag(-h/2*m3*ones(1,n))+diag(-h/2*m3*ones(1,n-1),1);
% mat4(2:2:end,2:6:end)=diag(-h/2*m1*ones(1,n))+diag(-h/2*m1*ones(1,n-1),1); mat4(2:2:end,4:6:end)=diag(-h/2*m2*ones(1,n))+diag(-h/2*m2*ones(1,n-1),1); mat4(2:2:end,6:6:end)=diag(-h/2*m3*ones(1,n))+diag(-h/2*m3*ones(1,n-1),1);
% mat4=mat4(1:(end-2),:);
% Aeq((6*(n-1)+1):(8*(n-1)),(21*n+1):(27*n))=mat3; Aeq((6*(n-1)+1):(8*(n-1)),(27*n+1):(33*n))=mat4;
% Aeq((8*(n-1)+1):(10*(n-1)),(27*n+1):(33*n))=mat3; Aeq((8*(n-1)+1):(10*(n-1)),(33*n+1):(39*n))=mat4;
% Collocation on P0
% mat3=zeros(2*n,8*n); mat3(1:2:(2*n-1),1:8:(8*n-7))=diag(-1*ones(1,n))+diag(ones(1,n-1),1); mat3(2:2:(2*n),2:8:(8*n-6))=diag(-1*ones(1,n))+diag(ones(1,n-1),1); mat3=mat3(1:(end-2),:); % -1,1 matrix, 2*(n-1) x 8*n
% mat4=diag(-1*ones(1,2*n))+diag(ones(1,2*n-2),2); mat4=mat4(1:(end-2),:); % -1,1 matrix, 2*(n-1) x 2*n
% mat5=diag(-h/2*ones(1,2*n))+diag(-h/2*ones(1,2*n-2),2); mat5=mat5(1:(end-2),:); %-h/2,-h/2 matrix, 2*(n-1) x 2*n
% Aeq((6*(n-1)+1):(8*(n-1)),(9*n+1):(17*n))=mat3; Aeq((6*(n-1)+1):(8*(n-1)),(17*n+1):(19*n))=mat5;
% Aeq((8*(n-1)+1):(10*(n-1)),(17*n+1):(19*n))=mat4; Aeq((8*(n-1)+1):(10*(n-1)),(19*n+1):(21*n))=mat5;
% Collocation on G3
% mat3=zeros(2*n,6*n); mat3(1:2:(2*n-1),5:6:6*n)=diag(-1*ones(1,n))+diag(ones(1,n-1),1); mat3(2:2:(2*n),6:6:6*n)=diag(-1*ones(1,n))+diag(ones(1,n-1),1); mat3=mat3(1:(end-2),:); % -1,1 matrix, 2*(n-1) x 6*n
% mat4=zeros(2*n,6*n); mat4(1:2:(2*n-1),5:6:6*n)=diag(-h/2*ones(1,n))+diag(-h/2*ones(1,n-1),1); mat4(2:2:(2*n),6:6:6*n)=diag(-h/2*ones(1,n))+diag(-h/2*ones(1,n-1),1); mat4=mat4(1:(end-2),:); % -h/2,-h/2 matrix, 2*(n-1) x 6*n
% Aeq((6*(n-1)+1):(8*(n-1)),(21*n+1):(27*n))=mat3; Aeq((6*(n-1)+1):(8*(n-1)),(27*n+1):(33*n))=mat4;
% Aeq((8*(n-1)+1):(10*(n-1)),(27*n+1):(33*n))=mat3; Aeq((8*(n-1)+1):(10*(n-1)),(33*n+1):(39*n))=mat4;
    %beq for collocation constraints are all 0.

% Boundary constraints
Aeq_bc=zeros(n_bc,n_dv); beq_bc=zeros(n_bc,1);
Aeq_bc(1:3,(3*n+1):(3*n+3))=diag(ones(1,3)); %beq=0
Aeq_bc(4:6,(6*n-2):(6*n))=diag(ones(1,3)); %beq=0
Aeq_bc(7:9,(6*n+1):(6*n+3))=diag(ones(1,3)); %beq=0
Aeq_bc(10:12,(9*n-2):(9*n))=diag(ones(1,3)); %beq=0
Aeq_bc(13:14,(9*n+1):(9*n+2))=diag(ones(1,2)); %beq=0
Aeq_bc(15:16,(17*n-7):(17*n-6))=diag(ones(1,2)); beq_bc(15:16)=[d_jump,h_jump];
Aeq_bc(17:19,1:3)=diag(ones(1,3)); beq_bc(17:19)=qi_list;
Aeq_bc(20:22,(3*n-2):(3*n))=diag(ones(1,3)); beq_bc(20:22)=qf_list;
Aeq_bc(23:24,(17*n+1):(17*n+2))=diag(ones(1,2)); %beq=0 %initial vP0ij=0
Aeq_bc(25:26,(19*n-1):(19*n))=diag(ones(1,2)); %beq=0 %final vP0ij=0
Aeq_bc(27:28,(19*n+1):(19*n+2))=diag(ones(1,2)); %beq=0 %initial aP0ij=0
Aeq_bc(29:30,(21*n-1):(21*n))=diag(ones(1,2)); %beq=0 %final aP0ij=0
Aeq((n_cl+1):(n_cl+n_bc),:)=Aeq_bc; beq((n_cl+1):(n_cl+n_bc))=beq_bc;

% Inequality constraints:
n_hu=4*n;
A=zeros(n_hu,n_dv); b=zeros(n_hu,1);
A(1:n,3:3:3*n)=diag(ones(1,n)); A(1:n,2:3:3*n)=diag(-1*ones(1,n)); b(1:n)=150/180*pi; %q3-q2<=150
A((n+1):2*n,3:3:3*n)=diag(-1*ones(1,n)); A((n+1):2*n,2:3:3*n)=diag(ones(1,n)); %b=0; %q3-q2>=0
A((2*n+1):3*n,1:3:3*n)=diag(ones(1,n)); A((2*n+1):3*n,2:3:3*n)=diag(-1*ones(1,n)); b((2*n+1):3*n)=150/180*pi; %q2-q1<=150
A((3*n+1):4*n,1:3:3*n)=diag(-1*ones(1,n)); A((3*n+1):4*n,2:3:3*n)=diag(ones(1,n)); %b=0; %q2-q1>=0

% ------------------------------------------------------------------------

% Initialization
var_list_guess=initialize_standing_jump_contact_implicit(n,qi_list,qf_list,L_list,d_list,m_list,I_list,g,d_jump,h_jump,tf);

% Optimization with fmincon
func_cost=@(var_list)standing_jump_cost(var_list,n,tf);
func_nlcon=@(var_list)standing_jump_nlcon_contact_implicit(var_list,n,m_list,I_list,g,L_list,d_list,d_jump,h_jump,miu,tf,h);
options = optimoptions('fmincon','Algorithm','interior-point','SubproblemAlgorithm','cg','MaxFunctionEvaluations',100e4,'MaxIterations',1000,'StepTol',1e-8,'PlotFcn',{'optimplotfval','optimplotconstrviolation'}); % 'interior-point' 'factorization'
[var_list,cost]=fmincon(func_cost,var_list_guess,A,b,Aeq,beq,[],[],func_nlcon,options);

toc;

fprintf('\nCost of the optimal path: %g.\n',cost);

% Visualize
P_list=var_list((9*n+1):(17*n)); 
P_mat=reshape(P_list,8,n);
figure(61);
for k=1:n
    plot(P_mat(1:2:end,k)',P_mat(2:2:end,k)','-bo'); hold on;
    plot([-0.5,d_jump-0.1,d_jump-0.1,d_jump+0.5],[0,0,h_jump,h_jump]); hold off;
    %plot([d_jump-0.5,d_jump,d_jump,0],[h_jump,h_jump,0,0]); hold off;
    axis([-0.5,d_jump+0.5,0,h_jump+1.6]);
    %axis([-0.5,d_jump+0.5,0,h_jump+5]);
    dt=tf/n; % Finish animation in tf
    pause(dt);
    if k==1
        pause(0.5);
    end
end

u_mat=reshape(var_list((39*n+1):(42*n)),3,n);
effort_list=sum(u_mat.^2);
% figure(62);
% plot(t_list,effort_list,'-r');
% xlabel('time(s)'); ylabel('Sum torque squared');

figure(63);
plot(t_list,u_mat(1,:)); hold on; plot(t_list,u_mat(2,:)); plot(t_list,u_mat(3,:)); hold off
xlabel('time(s)'); ylabel('Torque');
legend('u1','u2','u3');

clear all; clc;

tic;

% Decision variables:
% q(1,2,3,4,5),vq(1,2,3,4,5),aq(1,2,3,4,5),P(i,j)(1,2,3,4,5),G(i,j)(1,2,3,4,5),vG(i,j)(1,2,3,4,5),aG(i,j)(1,2,3,4,5),u(1,2,3,4,5)

n=10; % Number of knot points including 2 ends
t1=0; tf=1; t_list=linspace(t1,tf,n); h=(tf-t1)/(n-1); % Step size
qi_list=[-10,-20,0,10,20]/180*pi; qf_list=fliplr(qi_list);
d_step=0.3;
v_walk=d_step/tf;

% Model parameters:
g=9.81;
m_list=[3.2,6.8,20,6.8,3.2]; I_list=[0.93,1.08,2.22,1.08,0.93]; L_list=[0.4,0.4,0.625,0.4,0.4]; d_list=[0.272,0.237,0.2,0.163,0.128];

% Linear constraints------------------------------------------------------
n_cl=10*(n-1); % collocation between (q,vq) and (vq,aq)
%n_bc=7; % 5 for initial q, 2 for P5ij(n-1)=P5ij(n)
n_bc=4; % 2 for P5ij(n)=(D,0), 2 for P5ij(n-1)=P5ij(n)
n_rp=15; % repeatability: 2nd start = 1st start for all q,vq,aq but flipped in geometric order
Aeq=zeros(n_cl+n_bc+n_rp,59*n); beq=zeros(n_cl+n_bc+n_rp,1);

% Collocation constraints (trapezoidal)
mat1=diag(-1*ones(1,5*n))+diag(ones(1,5*(n-1)),5); mat1=mat1(1:(end-5),:); % -1,1 matrix, 5*(n-1) x 5*n
mat2=diag(-h/2*ones(1,5*n))+diag(-h/2*ones(1,5*(n-1)),5); mat2=mat2(1:(end-5),:); %-h/2,-h/2 matrix, 5*(n-1) x 5*n
Aeq(1:(5*(n-1)),1:(5*n))=mat1; Aeq(1:(5*(n-1)),(5*n+1):(10*n))=mat2; % q12345 % vq12345
Aeq((5*(n-1)+1):(10*(n-1)),(5*n+1):(10*n))=mat1; Aeq((5*(n-1)+1):(10*(n-1)),(10*n+1):(15*n))=mat2; % vq12345 & aq12345
    %beq for collocation constraints are all 0.
% mat1=diag(-1*ones(1,5*n))+diag(ones(1,5*(n-1)),5); mat1=mat1(1:(end-5),:); % -1,1 matrix, 5*(n-1) x 5*n
% mat2=diag(-h/2*ones(1,5*n))+diag(-h/2*ones(1,5*(n-1)),5); mat2=mat2(1:(end-5),:); %-h/2,-h/2 matrix, 5*(n-1) x 5*n
% mat11=diag(-1*ones(1,5*n))+diag(ones(1,5*(n-1)),5); mat11=mat11(1:(end-10),:); % -1,1 matrix, 5*(n-1)-1 x 5*n
% mat22=diag(-h/2*ones(1,5*n))+diag(-h/2*ones(1,5*(n-1)),5); mat22=mat22(1:(end-10),:); %-h/2,-h/2 matrix, 5*(n-1)-1 x 5*n
% Aeq(1:(5*(n-1)),1:(5*n))=mat1; Aeq(1:(5*(n-1)),(5*n+1):(10*n))=mat2; % q12345 % vq12345
% Aeq((5*(n-1)+1):(5*(n-1)+5*(n-2)),(5*n+1):(10*n))=mat11; Aeq((5*(n-1)+1):(5*(n-1)+5*(n-2)),(10*n+1):(15*n))=mat22; % vq12345 & aq12345 % exception between vaq5n and vaq5(n-1)
% mat1=diag(-1*ones(1,5*n))+diag(ones(1,5*(n-1)),5); mat1=mat1(1:(end-5),:); % -1,1 matrix, 5*(n-1) x 5*n
% mat2=diag(-h/2*ones(1,5*n))+diag(-h/2*ones(1,5*(n-1)),5); mat2=mat2(1:(end-5),:); %-h/2,-h/2 matrix, 5*(n-1) x 5*n
% mat11=diag(-1*ones(1,5*n))+diag(ones(1,5*(n-1)),5); mat11=mat11(1:(end-6),:); % -1,1 matrix, 5*(n-1)-1 x 5*n
% mat22=diag(-h/2*ones(1,5*n))+diag(-h/2*ones(1,5*(n-1)),5); mat22=mat22(1:(end-6),:); %-h/2,-h/2 matrix, 5*(n-1)-1 x 5*n
% Aeq(1:(5*(n-1)),1:(5*n))=mat1; Aeq(1:(5*(n-1)),(5*n+1):(10*n))=mat2; % q12345 % vq12345
% Aeq((5*(n-1)+1):(10*(n-1)-1),(5*n+1):(10*n))=mat11; Aeq((5*(n-1)+1):(10*(n-1)-1),(10*n+1):(15*n))=mat22; % vq12345 & aq12345 % exception between vaq5n and vaq5(n-1)


% Boundary constraints
% Aeq_bc=zeros(n_bc,59*n); beq_bc=zeros(n_bc,1);
% Aeq_bc(1:5,1:5)=eye(5); beq_bc(1:5)=qi_list';
% Aeq_bc(6,[(25*n-10),(25*n)])=[-1,1]; % beq=0;
% Aeq_bc(7,[(25*n-11),(25*n-1)])=[-1,1]; % beq=0;
% Aeq((n_cl+1):(n_cl+n_bc),:)=Aeq_bc; beq((n_cl+1):(n_cl+n_bc))=beq_bc;
Aeq_bc=zeros(n_bc,59*n); beq_bc=zeros(n_bc,1);
Aeq_bc(1,(25*n-1))=1; beq_bc(1)=d_step;
Aeq_bc(2,25*n)=1; beq_bc(2)=0;
Aeq_bc(3,[(25*n-10),(25*n)])=[-1,1]; % beq=0;
Aeq_bc(4,[(25*n-11),(25*n-1)])=[-1,1]; % beq=0;
Aeq((n_cl+1):(n_cl+n_bc),:)=Aeq_bc; beq((n_cl+1):(n_cl+n_bc))=beq_bc;

% Repeatability constraints: (av)q-(12345)=(av)q+(54321)
Aeq_rp=zeros(n_rp,59*n); beq_rp=zeros(n_rp,1);
Aeq_rp(1:5,1:5)=eye(5); Aeq_rp(1:5,(5*n-4):(5*n))=-fliplr(eye(5)); %beq=0
Aeq_rp(6:10,(5*n+1):(5*n+5))=eye(5); Aeq_rp(6:10,(10*n-4):(10*n))=-fliplr(eye(5)); %beq=0
Aeq_rp(11:15,(10*n+1):(10*n+5))=eye(5); Aeq_rp(11:15,(15*n-4):(15*n))=-fliplr(eye(5)); %beq=0
Aeq((n_cl+n_bc+1):(n_cl+n_bc+n_rp),:)=Aeq_rp; beq((n_cl+n_bc+1):(n_cl+n_bc+n_rp))=beq_rp;

% Inequality constraints:
n_path=n-3; n_human=4*n;
% Path constraints: P5j>=0.01 for path except end points
h_clearance=0.005;
A=zeros(n_path+n_human,59*n); b=zeros(n_path+n_human,1);
A(1:n_path,(15*n+20):10:(25*n-20))=-1*eye(n_path); b(1:n_path)=-h_clearance;

% Humanoid constraints
q_body=5/180*pi;
A_human=zeros(n_human,59*n);b_human=zeros(n_human,1);
A_human(1:n,3:5:(5*n))=eye(n); b_human(1:n)=q_body; %q3<q_body
A_human((n+1):(2*n),3:5:(5*n))=-eye(n); b_human((n+1):(2*n))=q_body; %-q3<q_body
A_human((2*n+1):(3*n),1:5:(5*n))=-eye(n); A_human((2*n+1):(3*n),2:5:(5*n))=eye(n); %b=0; %q2<q1
A_human((3*n+1):(4*n),4:5:(5*n))=eye(n); A_human((3*n+1):(4*n),5:5:(5*n))=-eye(n); %b=0; %q4<q5
A((n_path+1):(n_path+n_human),:)=A_human; b((n_path+1):(n_path+n_human))=b_human;

% ------------------------------------------------------------------------

% Initialization
var_list_guess=initialize_step_strike(n,t_list,qi_list,qf_list,L_list,d_list,m_list,I_list,g);

% Optimization with fmincon
options = optimoptions('fmincon','MaxFunctionEvaluations',20e4,'StepTol',1e-6);
func_cost=@(var_list)step_strike_cost(var_list,n,tf);
func_nlcon=@(var_list)step_strike_nlcon(var_list,n,m_list,I_list,g,L_list,d_list);
[var_list,cost]=fmincon(func_cost,var_list_guess,A,b,Aeq,beq,[],[],func_nlcon,options);

toc;

fprintf('\nCost of the optimal path: %g.\n',cost);

% Visualize
P_list=var_list((15*n+1):(25*n)); 
P_mat=reshape(P_list,10,n);
d_step=P_list(end-1);
n_step=10;
x_mat=zeros(7,n_step*(n-1)+1);
y_mat=zeros(7,n_step*(n-1)+1);
for i=1:n_step
    x_mat(:,((i-1)*(n-1)+1):(i*(n-1)))=[zeros(1,(n-1));P_mat([1,3,5,3,7,9],1:(n-1))]+(i-1)*d_step;
    y_mat(:,((i-1)*(n-1)+1):(i*(n-1)))=[zeros(1,(n-1));P_mat([2,4,6,4,8,10],1:(n-1))];
end
x_mat(:,end)=[0;P_mat([1,3,5,3,7,9],n)]+(n_step-1)*d_step;
y_mat(:,end)=[0;P_mat([2,4,6,4,8,10],n)];
figure(24);
for k=1:size(x_mat,2)
    plot(x_mat(:,k),y_mat(:,k),'-bo');
    axis([-1,1+(n_step-1)*d_step,0,1.5]);
    title('Periodic Walking Animation','FontSize', 15);
    xlabel('x (m)','FontSize', 15); ylabel('z (m)','FontSize', 15);
    dt=tf/n;
    pause(dt);
end

% figure(24);
% for i=1:n_step
%     for k=1:n
%         plot([0,P_mat([1,3,5,3,7,9],k)']+(i-1)*d_step,[0,P_mat([2,4,6,4,8,10],k)'],'-bo');
%         axis([-1,1+(n_step-1)*d_step,0,1.6]);
%         dt=tf/n; % Finish animation in tf
%         pause(dt);
%     end
% end

u_mat=reshape(var_list((55*n+1):(59*n)),4,n);
effort_list=sum(u_mat.^2);
figure(25);
plot(t_list,effort_list,'-r');
xlabel('time(s)'); ylabel('Sum torque squared');

figure(26);
plot(t_list,u_mat(1,:)); hold on; plot(t_list,u_mat(2,:)); plot(t_list,u_mat(3,:)); plot(t_list,u_mat(4,:)); hold off
xlabel('time(s)'); ylabel('Torque');
legend('u2','u3','u4','u5');

P_list=var_list((15*n+1):(25*n)); 
P_mat=reshape(P_list,10,n);
figure(27);
plot(P_mat(9,:),P_mat(10,:));
xlabel('P5i'); ylabel('P5j');



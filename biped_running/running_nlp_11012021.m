clear all; clc;

tic;

t1=0; tf=1; d_run=0.4; h_run=0.3; iswalking=true;

n=15; % Number of knot points including 2 ends
% Decision variables:
% q(5*n),vq(5*n),aq(5*n),P(12*n),vP0(2*n),aP0(2*n),vP5(2*n),aP5(2*n),G(10*n),vG(10*n),aG(10*n),GRF(4*n),u(6*n),t1,t2
n_dv=75*n+2; % (+2 time variables for phase transition)

t_list=linspace(t1,tf,n); h=(tf-t1)/(n-1); % Step size
qi_list=[5,-10,5,5,10]/180*pi; qf_list=fliplr(qi_list);
miu=0.8; % coefficient of friction
ft_clr=0;

% Model parameters:
g=9.81;
m_list=[3.2,6.8,20,6.8,3.2]; I_list=[0.93,1.08,2.22,1.08,0.93]; L_list=[0.4,0.4,0.625,0.4,0.4]; d_list=[0.128,0.163,0.2,0.163,0.128];

% Linear constraints------------------------------------------------------
n_cl=10*(n-1); % collocation between (q,vq) and (vq,aq)
n_bc=13; % 2 for P5ij(n)=(D,0), 2 for P5ij(n-1)=P5ij(n)
n_rp=21;% repeatability
n_GRF=2*n;
Aeq=zeros(n_cl+n_bc+n_rp+n_GRF,n_dv); beq=zeros(n_cl+n_bc+n_rp+n_GRF,1);

% Collocation constraints (trapezoidal)
mat1=diag(-1*ones(1,5*n))+diag(ones(1,5*(n-1)),5); mat1=mat1(1:(end-5),:); % -1,1 matrix, 5*(n-1) x 5*n
mat2=diag(-h/2*ones(1,5*n))+diag(-h/2*ones(1,5*(n-1)),5); mat2=mat2(1:(end-5),:); %-h/2,-h/2 matrix, 5*(n-1) x 5*n
Aeq(1:(5*(n-1)),1:(5*n))=mat1; Aeq(1:(5*(n-1)),(5*n+1):(10*n))=mat2; % q12345 % vq12345
Aeq((5*(n-1)+1):(10*(n-1)),(5*n+1):(10*n))=mat1; Aeq((5*(n-1)+1):(10*(n-1)),(10*n+1):(15*n))=mat2; % vq12345 & aq12345
    %beq for collocation constraints are all 0.

% Boundary constraints
Aeq_bc=zeros(n_bc,n_dv); beq_bc=zeros(n_bc,1);
Aeq_bc(1:2,(15*n+1):(15*n+2))=eye(2); %beq_bc=0; % initial P0ij=0
Aeq_bc(3:4,(27*n+1):(27*n+2))=eye(2); %beq_bc=0; % initial vP0ij=0
Aeq_bc(5:6,(29*n+1):(29*n+2))=eye(2); %beq_bc=0; % initial aP0ij=0
Aeq_bc(7:8,(27*n-1):(27*n))=eye(2); beq_bc(7:8)=[d_run,h_run]; % final P5ij
Aeq_bc(9:10,(33*n-1):(33*n))=eye(2); %beq_bc=0; % final vP5ij=0
Aeq_bc(11:12,(35*n-1):(35*n))=eye(2); %beq_bc=0; % final aP5ij=0
if iswalking
    Aeq_bc(13,75*n+1)=1; Aeq_bc(13,75*n+2)=-1; %beq=0; % t1=t2 for walking 
end
Aeq((n_cl+1):(n_cl+n_bc),:)=Aeq_bc; beq((n_cl+1):(n_cl+n_bc))=beq_bc;

% Repeatability constraints: (av)q-(12345)=(av)q+(54321)
Aeq_rp=zeros(n_rp,n_dv); beq_rp=zeros(n_rp,1);
Aeq_rp(1:5,1:5)=eye(5); Aeq_rp(1:5,(5*n-4):(5*n))=-fliplr(eye(5)); %beq=0
Aeq_rp(6:10,(5*n+1):(5*n+5))=eye(5); Aeq_rp(6:10,(10*n-4):(10*n))=-fliplr(eye(5)); %beq=0
Aeq_rp(11:15,(10*n+1):(10*n+5))=eye(5); Aeq_rp(11:15,(15*n-4):(15*n))=-fliplr(eye(5)); %beq=0
Aeq_rp(16:21,(69*n+1):(69*n+6))=eye(6); Aeq_rp(16:21,(75*n-5):(75*n))=fliplr(eye(6)); %beq=0  % initial u123456 = - final u654321

Aeq((n_cl+n_bc+1):(n_cl+n_bc+n_rp),:)=Aeq_rp; beq((n_cl+n_bc+1):(n_cl+n_bc+n_rp))=beq_rp;

% !!! add repeatability of u !!!

% Ground Reaction Force constraints: GRF0ij+GRF5ij-sum(m*g*j)-sum(m*aGij)=0
Aeq_GRF=zeros(n_GRF,n_dv); beq_GRF=zeros(n_GRF,1);
mat_m=zeros(n,5*n);
for i=1:5
    mat_m(1:n,i:5:5*n)=m_list(i)*eye(n);
end
Aeq_GRF(1:n,(65*n+1):4:(69*n))=eye(n); Aeq_GRF(1:n,(65*n+3):4:(69*n))=eye(n); Aeq_GRF(1:n,(55*n+1):2:(65*n))=-mat_m; %beq_GRF=0; % i
Aeq_GRF((n+1):2*n,(65*n+2):4:(69*n))=eye(n); Aeq_GRF((n+1):2*n,(65*n+4):4:(69*n))=eye(n); Aeq_GRF((n+1):2*n,(55*n+2):2:(65*n))=-mat_m; beq_GRF((n+1):2*n)=sum(m_list)*g; % j
Aeq((n_cl+n_bc+n_rp+1):(n_cl+n_bc+n_rp+n_GRF),:)=Aeq_GRF; beq((n_cl+n_bc+n_rp+1):(n_cl+n_bc+n_rp+n_GRF))=beq_GRF;

% Inequality constraints:
n_hu=4*n; % humanoid constraints
n_ph=3; % phase constraints
n_fN=6*n; % Ground reaction force inequality constraints
A=zeros(n_hu+n_ph+n_fN,n_dv); b=zeros(n_hu+n_ph+n_fN,1);
% Humanoid constraints
q_body=5/180*pi;
A(1:n,3:5:(5*n))=eye(n); b(1:n)=q_body; %q3<=q_body
A((n+1):(2*n),3:5:(5*n))=-eye(n); b((n+1):(2*n))=q_body; %-q3<=q_body
A((2*n+1):(3*n),1:5:(5*n))=-eye(n); A((2*n+1):(3*n),2:5:(5*n))=eye(n); %b=0; %q2<=q1
A((3*n+1):(4*n),4:5:(5*n))=eye(n); A((3*n+1):(4*n),5:5:(5*n))=-eye(n); %b=0; %q4<=q5
% Phase constraints
A_ph=zeros(n_ph,n_dv); b_ph=zeros(n_ph,1);
A_ph(1,75*n+1)=-1; b_ph(1)=-1e-2; % t1>0
A_ph(2,75*n+2)=1; b_ph(2)=tf-1e-2; % t2<tf
A_ph(3,75*n+1)=1; A_ph(3,75*n+2)=-1; % b=0; % t1<=t2
A((n_hu+1):(n_hu+n_ph),:)=A_ph; b((n_hu+1):(n_hu+n_ph))=b_ph;
% Ground reaction force inequality constraints
A_fN=zeros(n_fN,n_dv); b_fN=zeros(n_fN,1);
A_fN(1:n,(65*n+1):4:(69*n))=eye(n); A_fN(1:n,(65*n+2):4:(69*n))=-miu*eye(n); %b_fN=0 % f0<=miu*N0
A_fN((n+1):2*n,(65*n+1):4:(69*n))=-eye(n); A_fN((n+1):2*n,(65*n+2):4:(69*n))=-miu*eye(n); %b_fN=0 % -f0<=miu*N0
A_fN((2*n+1):3*n,(65*n+3):4:(69*n))=eye(n); A_fN((2*n+1):3*n,(65*n+4):4:(69*n))=-miu*eye(n); %b_fN=0 % f5<=miu*N5
A_fN((3*n+1):4*n,(65*n+3):4:(69*n))=-eye(n); A_fN((3*n+1):4*n,(65*n+4):4:(69*n))=-miu*eye(n); %b_fN=0 % -f5<=miu*N5
A_fN((4*n+1):5*n,(65*n+2):4:(69*n))=-eye(n); %b_fN=0 % N0>=0 (-N0<=0)
A_fN((5*n+1):6*n,(65*n+4):4:(69*n))=-eye(n); %b_fN=0 % N5>=0 (-N5<=0)
A((n_hu+n_ph+1):(n_hu+n_ph+n_fN),:)=A_fN; b((n_hu+n_ph+1):(n_hu+n_ph+n_fN))=b_fN;

% ------------------------------------------------------------------------

% Initialization
var_list_guess=initialize_running(n,qi_list,qf_list,L_list,d_list,m_list,I_list,g,d_run,h_run,t_list,iswalking);

% Optimization with fmincon
func_cost=@(var_list)running_cost(var_list,n,tf);
func_nlcon=@(var_list)running_nlcon(var_list,n,m_list,I_list,g,L_list,d_list,d_run,h_run,tf,h,ft_clr);
options = optimoptions('fmincon','Algorithm','interior-point','SubproblemAlgorithm','cg','MaxFunctionEvaluations',100e4,'MaxIterations',1000,'StepTol',1e-8,'PlotFcn',{'optimplotfval','optimplotconstrviolation'}); % 'interior-point' 'factorization'/'cg'
[var_list,cost]=fmincon(func_cost,var_list_guess,A,b,Aeq,beq,[],[],func_nlcon,options);

toc;

fprintf('\nCost of the optimal path: %g.\n',cost);

% Animation
P_list=var_list((15*n+1):(27*n)); 
P_mat=reshape(P_list,12,n);
n_step=5;
x_mat=zeros(7,n_step*(n-1));
y_mat=zeros(7,n_step*(n-1));
for i=1:n_step
    x_mat(:,((i-1)*(n-1)+1):(i*(n-1)))=[P_mat([1,3,5,7,5,9,11],1:(n-1))]+(i-1)*d_run;
    y_mat(:,((i-1)*(n-1)+1):(i*(n-1)))=[P_mat([2,4,6,8,6,10,12],1:(n-1))]+(i-1)*h_run;
end
figure(48);
for k=1:size(x_mat,2)
    plot(x_mat(:,k),y_mat(:,k),'-bo'); hold on;
    for kk=1:n_step
    plot([kk*d_run-0.1,kk*d_run-0.1,(kk+1)*d_run],[(kk-1)*h_run,kk*h_run,kk*h_run]); 
    end
    hold off;
    axis([-0.5,(n_step+1)*d_run,0,1.5+(n_step)*h_run]); set(gca,'DataAspectRatio',[1 1 1]);
    title('Running Animation','FontSize', 15);
    xlabel('x (m)','FontSize', 15); ylabel('z (m)','FontSize', 15);
    dt=tf/n/1000;
    pause(dt);
    if k==1
        pause(1);
    end
end

u_mat=reshape(var_list((69*n+1):(75*n)),6,n);
effort_list=sum(u_mat.^2);
% figure(62);
% plot(t_list,effort_list,'-r');
% xlabel('time(s)'); ylabel('Sum torque squared');

figure(58);
plot(t_list,u_mat(1,:)); hold on; plot(t_list,u_mat(2,:)); plot(t_list,u_mat(3,:)); plot(t_list,u_mat(4,:)); plot(t_list,u_mat(5,:)); plot(t_list,u_mat(6,:)); hold off;
xlabel('time(s)'); ylabel('Torque');
legend('u1','u2','u3','u4','u5','u6');

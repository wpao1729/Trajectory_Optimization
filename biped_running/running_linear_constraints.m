function [Aeq,beq,A,b] = running_linear_constraints(n,m_list,g,d_run,h_run,tf,h,n_dv,iswalking,miu)

% Linear equality constraints------------------------------------------------------
n_cl=10*(n-1); % number of collocation constraints between (q,vq) and (vq,aq)
n_bc=2; % number of boundary constraints
n_rp=21;% repeatability
n_GRF=2*n;
Aeq=zeros(n_cl+n_bc+n_rp+n_GRF,n_dv); beq=zeros(n_cl+n_bc+n_rp+n_GRF,1);

% Collocation constraints (trapezoidal collocation)
mat1=diag(-1*ones(1,5*n))+diag(ones(1,5*(n-1)),5); mat1=mat1(1:(end-5),:); % -1,1 matrix, 5*(n-1) x 5*n
mat2=diag(-h/2*ones(1,5*n))+diag(-h/2*ones(1,5*(n-1)),5); mat2=mat2(1:(end-5),:); %-h/2,-h/2 matrix, 5*(n-1) x 5*n
Aeq(1:(5*(n-1)),1:(5*n))=mat1; Aeq(1:(5*(n-1)),(5*n+1):(10*n))=mat2; % collocation between q12345 & vq12345
Aeq((5*(n-1)+1):(10*(n-1)),(5*n+1):(10*n))=mat1; Aeq((5*(n-1)+1):(10*(n-1)),(10*n+1):(15*n))=mat2; % collocation between vq12345 & aq12345
    %beq for collocation constraints are all 0.

% Boundary constraints (t1 and t2)
Aeq_bc=zeros(n_bc,n_dv); beq_bc=zeros(n_bc,1);
if iswalking
    Aeq_bc(1,75*n+1)=1; Aeq_bc(1,75*n+2)=-1; %beq=0; % t1=t2 for walking 
end
Aeq_bc(2,75*n+1)=1; Aeq_bc(2,75*n+2)=1; beq_bc(2)=tf; % t1+t2=tf
Aeq((n_cl+1):(n_cl+n_bc),:)=Aeq_bc; beq((n_cl+1):(n_cl+n_bc))=beq_bc;

% Repeatability constraints: (av)q-(12345)=(av)q+(54321)
Aeq_rp=zeros(n_rp,n_dv); beq_rp=zeros(n_rp,1);
Aeq_rp(1:5,1:5)=eye(5); Aeq_rp(1:5,(5*n-4):(5*n))=-fliplr(eye(5)); %beq=0 % initial q12345 = final q54321
Aeq_rp(6:10,(5*n+1):(5*n+5))=eye(5); Aeq_rp(6:10,(10*n-4):(10*n))=-fliplr(eye(5)); %beq=0 % initial vq12345 = final vq54321
Aeq_rp(11:15,(10*n+1):(10*n+5))=eye(5); Aeq_rp(11:15,(15*n-4):(15*n))=-fliplr(eye(5)); %beq=0 % initial aq12345 = final aq54321
Aeq_rp(16:21,(69*n+1):(69*n+6))=eye(6); Aeq_rp(16:21,(75*n-5):(75*n))=fliplr(eye(6)); %beq=0  % initial u123456 = - final u654321

Aeq((n_cl+n_bc+1):(n_cl+n_bc+n_rp),:)=Aeq_rp; beq((n_cl+n_bc+1):(n_cl+n_bc+n_rp))=beq_rp;

% Ground Reaction Force constraints: GRF0ij+GRF5ij-sum(m*g*j)-sum(m*aGij)=0
Aeq_GRF=zeros(n_GRF,n_dv); beq_GRF=zeros(n_GRF,1);
mat_m=zeros(n,5*n);
for i=1:5
    mat_m(1:n,i:5:5*n)=m_list(i)*eye(n);
end
Aeq_GRF(1:n,(65*n+1):4:(69*n))=eye(n); Aeq_GRF(1:n,(65*n+3):4:(69*n))=eye(n); Aeq_GRF(1:n,(55*n+1):2:(65*n))=-mat_m; %beq_GRF=0; % GRF0i + GRF5i = sum(m*aGi)
Aeq_GRF((n+1):2*n,(65*n+2):4:(69*n))=eye(n); Aeq_GRF((n+1):2*n,(65*n+4):4:(69*n))=eye(n); Aeq_GRF((n+1):2*n,(55*n+2):2:(65*n))=-mat_m; beq_GRF((n+1):2*n)=sum(m_list)*g; % GRF0j + GRF5j - sum(m*gj) = sum(m*aGj)
Aeq((n_cl+n_bc+n_rp+1):(n_cl+n_bc+n_rp+n_GRF),:)=Aeq_GRF; beq((n_cl+n_bc+n_rp+1):(n_cl+n_bc+n_rp+n_GRF))=beq_GRF;

% Linear inequality constraints --------------------------------------------------
n_hu=4*n; % number of humanoid constraints
n_ph=3; % number of phase constraints
n_fN=6*n; % number of Ground Reaction Force inequality constraints
A=zeros(n_hu+n_ph+n_fN,n_dv); b=zeros(n_hu+n_ph+n_fN,1);

% Humanoid constraints
% Limit upper body position to -5 ~ 5 degree
q_body=5/180*pi; 
A(1:n,3:5:(5*n))=eye(n); b(1:n)=q_body; % q3<=q_body
A((n+1):(2*n),3:5:(5*n))=-eye(n); b((n+1):(2*n))=q_body; % -q3<=q_body
% Knee cannot bend backward (over 180 degree)
A((2*n+1):(3*n),1:5:(5*n))=-eye(n); A((2*n+1):(3*n),2:5:(5*n))=eye(n); %b=0; %q2<=q1
A((3*n+1):(4*n),4:5:(5*n))=eye(n); A((3*n+1):(4*n),5:5:(5*n))=-eye(n); %b=0; %q4<=q5

% Phase constraints (t1 and t2)
A_ph=zeros(n_ph,n_dv); b_ph=zeros(n_ph,1);
A_ph(1,75*n+1)=-1; b_ph(1)=-1e-2; % t1>0
A_ph(2,75*n+2)=1; b_ph(2)=tf-1e-2; % t2<tf
A_ph(3,75*n+1)=1; A_ph(3,75*n+2)=-1; % b=0; % t1<=t2
A((n_hu+1):(n_hu+n_ph),:)=A_ph; b((n_hu+1):(n_hu+n_ph))=b_ph;

% Ground reaction force inequality constraints
A_fN=zeros(n_fN,n_dv); b_fN=zeros(n_fN,1);
% static friction constraint: |f|<=miu*N
A_fN(1:n,(65*n+1):4:(69*n))=eye(n); A_fN(1:n,(65*n+2):4:(69*n))=-miu*eye(n); %b_fN=0 % f0<=miu*N0
A_fN((n+1):2*n,(65*n+1):4:(69*n))=-eye(n); A_fN((n+1):2*n,(65*n+2):4:(69*n))=-miu*eye(n); %b_fN=0 % -f0<=miu*N0
A_fN((2*n+1):3*n,(65*n+3):4:(69*n))=eye(n); A_fN((2*n+1):3*n,(65*n+4):4:(69*n))=-miu*eye(n); %b_fN=0 % f5<=miu*N5
A_fN((3*n+1):4*n,(65*n+3):4:(69*n))=-eye(n); A_fN((3*n+1):4*n,(65*n+4):4:(69*n))=-miu*eye(n); %b_fN=0 % -f5<=miu*N5
% Normal ground reaction force must be >=0
A_fN((4*n+1):5*n,(65*n+2):4:(69*n))=-eye(n); %b_fN=0 % N0>=0 (-N0<=0)
A_fN((5*n+1):6*n,(65*n+4):4:(69*n))=-eye(n); %b_fN=0 % N5>=0 (-N5<=0)
A((n_hu+n_ph+1):(n_hu+n_ph+n_fN),:)=A_fN; b((n_hu+n_ph+1):(n_hu+n_ph+n_fN))=b_fN;


end
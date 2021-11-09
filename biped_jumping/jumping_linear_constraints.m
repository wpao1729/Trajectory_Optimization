function [Aeq,beq,A,b] = jumping_linear_constraints (n,d_jump,h_jump,tf,h,qi_list,qf_list,n_dv)

% Linear equality constraints----------------------------------------------
n_cl=6*(n-1); % number of collocation constraints between (q,vq) and (vq,aq)
n_bc=22; % number of boundary constraints
Aeq=zeros(n_cl+n_bc,n_dv); beq=zeros(n_cl+n_bc,1);

% Collocation constraints (trapezoidal)
mat1=diag(-1*ones(1,3*n))+diag(ones(1,3*(n-1)),3); mat1=mat1(1:(end-3),:); % -1,1 matrix, 3*(n-1) x 3*n
mat2=diag(-h/2*ones(1,3*n))+diag(-h/2*ones(1,3*(n-1)),3); mat2=mat2(1:(end-3),:); %-h/2,-h/2 matrix, 3*(n-1) x 3*n
Aeq(1:(3*(n-1)),1:(3*n))=mat1; Aeq(1:(3*(n-1)),(3*n+1):(6*n))=mat2; % collocation between q123 % vq123
Aeq((3*(n-1)+1):(6*(n-1)),(3*n+1):(6*n))=mat1; Aeq((3*(n-1)+1):(6*(n-1)),(6*n+1):(9*n))=mat2; % collocation between vq123 & aq123
    %beq for collocation constraints are all 0.

% Boundary constraints
Aeq_bc=zeros(n_bc,n_dv); beq_bc=zeros(n_bc,1);
Aeq_bc(1:3,(3*n+1):(3*n+3))=diag(ones(1,3)); %beq=0 % initial vq123 = 0
Aeq_bc(4:6,(6*n-2):(6*n))=diag(ones(1,3)); %beq=0 % final vq123 = 0
Aeq_bc(7:9,(6*n+1):(6*n+3))=diag(ones(1,3)); %beq=0 % initial aq123 = 0
Aeq_bc(10:12,(9*n-2):(9*n))=diag(ones(1,3)); %beq=0 % final aq123 = 0
Aeq_bc(13:14,(9*n+1):(9*n+2))=diag(ones(1,2)); %beq=0 % initial P0ij
Aeq_bc(15:16,(17*n-7):(17*n-6))=diag(ones(1,2)); beq_bc(15:16)=[d_jump,h_jump]; % final P0ij
Aeq_bc(17:19,1:3)=diag(ones(1,3)); beq_bc(17:19)=qi_list; % initial q123
Aeq_bc(20:22,(3*n-2):(3*n))=diag(ones(1,3)); beq_bc(20:22)=qf_list; % final q123
Aeq((n_cl+1):(n_cl+n_bc),:)=Aeq_bc; beq((n_cl+1):(n_cl+n_bc))=beq_bc;

% Linear inequality constraints -------------------------------------------
n_hu=6*n; % number of humanoid constraints
n_ph=3; % number of phase constraints (t1,t2)
A=zeros(n_hu+n_ph,n_dv); b=zeros(n_hu+n_ph,1);

% Humanoid constraints
A(1:n,3:3:3*n)=diag(ones(1,n)); A(1:n,2:3:3*n)=diag(-1*ones(1,n)); b(1:n)=150/180*pi; %q3-q2<=150
A((n+1):2*n,3:3:3*n)=diag(-1*ones(1,n)); A((n+1):2*n,2:3:3*n)=diag(ones(1,n)); %b=0; %q3-q2>=0
A((2*n+1):3*n,1:3:3*n)=diag(ones(1,n)); A((2*n+1):3*n,2:3:3*n)=diag(-1*ones(1,n)); b((2*n+1):3*n)=150/180*pi; %q2-q1<=150
A((3*n+1):4*n,1:3:3*n)=diag(-1*ones(1,n)); A((3*n+1):4*n,2:3:3*n)=diag(ones(1,n)); %b=0; %q2-q1>=0

% Phase constraints
A_ph=zeros(n_ph,n_dv); b_ph=zeros(n_ph,1);
A_ph(1,42*n+1)=-1; b_ph(1)=-1e-6; % t1>0
A_ph(2,42*n+2)=1; b_ph(2)=tf-1e-6; % t2<tf
A_ph(3,42*n+1)=1; A_ph(3,42*n+2)=-1; b_ph(3)=-1e-6; % t1<t2
A((n_hu+1):(n_hu+n_ph),:)=A_ph; b((n_hu+1):(n_hu+n_ph))=b_ph;

end
function [Aeq,beq,A,b] = block_move_linear_constraints (n,n_dv,d,m,h)

% Linear equality constraints ---------------------------------------------
n_bc=4; % number of boundary constraints
n_cl=2*(n-1); % number of collocation constraints
n_dy=n; % number of dynamic constraints
Aeq=zeros(n_bc+n_cl+n_dy,n_dv); beq=zeros(n_bc+n_cl+n_dy,1);

% Boundary constraints
x1=0; xf=d; v1=0; vf=0;
Aeq(1,1)=1; Aeq(2,n)=1; Aeq(3,n+1)=1; Aeq(4,2*n)=1;
beq(1:4)=[x1,xf,v1,vf];

% Collocation constraints - trapezoidal rule: x(k+1)-x(k)=h/2*(u(k+1)+u(k))
mat1=diag(-1*ones(1,n))+diag(ones(1,n-1),1); mat1=mat1(1:(end-1),:);
mat2=diag(-h/2*ones(1,n))+diag(-h/2*ones(1,n-1),1); mat2=mat2(1:(end-1),:);
Aeq((n_bc+1):(n_bc+(n-1)),1:n)=mat1; Aeq((n_bc+1):(n_bc+(n-1)),(n+1):(2*n))=mat2; % collocation between x and v
Aeq((n_bc+n):(n_bc+n_cl),(n+1):2*n)=mat1; Aeq((n_bc+n):(n_bc+n_cl),(2*n+1):(3*n))=mat2; % collocation between v and a
    % beq for collocation constraints are all 0.

% Dynamics Constraints: u=ma
Aeq((n_bc+n_cl+1):(n_bc+n_cl+n_dy),(2*n+1):(3*n))=-m*eye(n); Aeq((n_bc+n_cl+1):(n_bc+n_cl+n_dy),(3*n+1):(4*n))=eye(n); % beq=0; % u-ma=0

% Linear inequality constraints ------------------------------------------
A=[]; b=[];

end
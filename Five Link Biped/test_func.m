% syms a va [2 1]
% x1=sin(a(1))+cos(a(2));
% x2=a(1)^2+3*a(2);
% x=[x1,x2]';
% vx=jacobian(x,a)*va;

syms r [2 1]
x=2*r;
subs(x,'r1',1)

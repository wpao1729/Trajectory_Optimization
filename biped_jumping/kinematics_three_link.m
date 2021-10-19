clear all; clc;

syms q vq aq L d [3 1]
syms P0i P0j vP0i vP0j aP0i aP0j

P1i=P0i+L1*sin(q1); P1j=P0j+L1*cos(q1);
P2i=P1i+L2*sin(q2); P2j=P1j+L2*cos(q2);
P3i=P2i+L3*sin(q3); P3j=P2j+L3*cos(q3);
P=[P0i,P0j,P1i,P1j,P2i,P2j,P3i,P3j]';

G1i=P0i+d1*sin(q1); G1j=P0j+d1*cos(q1);
G2i=P1i+d2*sin(q2); G2j=P1j+d2*cos(q2);
G3i=P2i+d3*sin(q3); G3j=P2j+d3*cos(q3);
G=[G1i,G1j,G2i,G2j,G3i,G3j]';

% vP=jacobian(P,q)*vq;
% aP=jacobian(vP,q)*vq+jacobian(vP,vq)*aq;
vG=jacobian(G,q)*vq+jacobian(G,P0i)*vP0i+jacobian(G,P0j)*vP0j;
aG=jacobian(vG,q)*vq+jacobian(vG,vq)*aq+jacobian(vG,P0i)*vP0i+jacobian(vG,P0j)*vP0j+jacobian(vG,vP0i)*aP0i+jacobian(vG,vP0j)*aP0j;

% save('kinematics.mat','P','vP','aP','G','vG','aG');


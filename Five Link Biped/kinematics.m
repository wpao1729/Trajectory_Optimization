clear all; clc;

% L=[0.5,0.5,0.5,0.5,0.5];
% d=[0.25,0.25,0.25,0.25,0.25];
% L1=L(1); L2=L(2); L3=L(3); L4=L(4); L5=L(5);
% d1=d(1); d2=d(2); d3=d(3); d4=d(4); d5=d(5);

syms q vq aq L d [5 1]

P1i=L1*sin(q1); P1j=L1*cos(q1);
P2i=P1i+L2*sin(q2); P2j=P1j+L2*cos(q2);
P3i=P2i+L3*sin(q3); P3j=P2j+L3*cos(q3);
P4i=P2i-L4*sin(q4); P4j=P2j-L4*cos(q4);
P5i=P4i-L5*sin(q5); P5j=P4j-L5*cos(q5);
P=[P1i,P1j,P2i,P2j,P3i,P3j,P4i,P4j,P5i,P5j]';

G1i=d1*sin(q1); G1j=d1*cos(q1);
G2i=P1i+d2*sin(q2); G2j=P1j+d2*cos(q2);
G3i=P2i+d3*sin(q3); G3j=P2j+d3*cos(q3);
G4i=P2i-d4*sin(q4); G4j=P2j-d4*cos(q4);
G5i=P4i-d5*sin(q5); G5j=P4j-d5*cos(q5);
G=[G1i,G1j,G2i,G2j,G3i,G3j,G4i,G4j,G5i,G5j]';

vP=jacobian(P,q)*vq;
aP=jacobian(vP,q)*vq+jacobian(vP,vq)*aq;
vG=jacobian(G,q)*vq;
aG=jacobian(vG,q)*vq+jacobian(vG,vq)*aq;
% save('kinematics.mat','P','vP','aP','G','vG','aG');


clear all; clc;

% L=[0.5,0.5,0.5,0.5,0.5];
% d=[0.25,0.25,0.25,0.25,0.25];
% L1=L(1); L2=L(2); L3=L(3); L4=L(4); L5=L(5);
% d1=d(1); d2=d(2); d3=d(3); d4=d(4); d5=d(5);

syms q vq aq L d [5 1];
syms P5i P5j;

P4i=P5i+L5*sin(q5); P4j=P5j+L5*cos(q5);
P2i=P4i+L4*sin(q4); P2j=P4j+L4*cos(q4);
P3i=P2i+L3*sin(q3); P3j=P2j+L3*cos(q3);
P1i=P2i-L2*sin(q2); P1j=P2j-L2*cos(q2);

G5i=P5i+(L5-d5)*sin(q5); G5j=P5j+(L5-d5)*cos(q5);
G4i=P4i+(L4-d4)*sin(q4); G4j=P4j+(L4-d4)*cos(q4);
G3i=P2i+d3*sin(q3); G3j=P2j+d3*cos(q3);
G2i=P2i-(L2-d2)*sin(q2); G2j=P2j-(L2-d2)*cos(q2);
G1i=P1i-(L1-d1)*sin(q1); G1j=P1j-(L1-d1)*cos(q1);

G=[G1i,G1j,G2i,G2j,G3i,G3j,G4i,G4j,G5i,G5j]';

vG=jacobian(G,q)*vq;
aG=jacobian(vG,q)*vq+jacobian(vG,vq)*aq;
% save('kinematics.mat','P','vP','aP','G','vG','aG');


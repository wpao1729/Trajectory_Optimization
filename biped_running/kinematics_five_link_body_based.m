clear all; clc;

syms q vq aq L d [5 1]
syms G3i G3j vG3i vG3j aG3i aG3j

P3i=G3i+(L3-d3)*sin(q3); P3j=G3j+(L3-d3)*cos(q3);
P2i=G3i-d3*sin(q3); P2j=G3j-d3*cos(q3);
P1i=P2i-L2*sin(q2); P1j=P2j-L2*cos(q2);
P0i=P1i-L1*sin(q1); P0j=P1j-L1*cos(q1);
P4i=P2i-L4*sin(q4); P4j=P2j-L4*cos(q4);
P5i=P4i-L5*sin(q5); P5j=P4j-L5*cos(q5);
P=[P0i,P0j,P1i,P1j,P2i,P2j,P3i,P3j,P4i,P4j,P5i,P5j]';

G2i=P2i-d2*sin(q2); G2j=P2j-d2*cos(q2);
G1i=P1i-d1*sin(q1); G1j=P1j-d1*cos(q1);
G4i=P2i-d4*sin(q4); G4j=P2j-d4*cos(q4);
G5i=P4i-d5*sin(q5); G5j=P4j-d5*cos(q5);
G=[G1i,G1j,G2i,G2j,G3i,G3j,G4i,G4j,G5i,G5j]';

vP=jacobian(P,q)*vq+jacobian(P,G3i)*vG3i+jacobian(P,G3j)*vG3j;
aP=jacobian(vP,q)*vq+jacobian(vP,vq)*aq+jacobian(vP,G3i)*vG3i+jacobian(vP,G3j)*vG3j+jacobian(vP,vG3i)*aG3i+jacobian(vP,vG3j)*aG3j;

vG=jacobian(G,q)*vq+jacobian(G,G3i)*vG3i+jacobian(G,G3j)*vG3j;
aG=jacobian(vG,q)*vq+jacobian(vG,vq)*aq+jacobian(vG,G3i)*vG3i+jacobian(vG,G3j)*vG3j+jacobian(vG,vG3i)*aG3i+jacobian(vG,vG3j)*aG3j;


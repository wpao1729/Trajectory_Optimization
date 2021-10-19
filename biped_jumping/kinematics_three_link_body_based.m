clear all; clc;

syms q vq aq L d [3 1]
syms G3i G3j vG3i vG3j aG3i aG3j

P3i=G3i+(L3-d3)*sin(q3); P3j=G3j+(L3-d3)*cos(q3);
P2i=G3i-d3*sin(q3); P2j=G3j-d3*cos(q3);
P1i=P2i-L2*sin(q2); P1j=P2j-L2*cos(q2);
P0i=P1i-L1*sin(q1); P0j=P1j-L1*cos(q1);
P=[P0i,P0j,P1i,P1j,P2i,P2j,P3i,P3j]';

G2i=P2i-(L2-d2)*sin(q2); G2j=P2j-(L2-d2)*cos(q2);
G1i=P1i-(L1-d1)*sin(q1); G1j=P1j-(L1-d1)*cos(q1);
G=[G1i,G1j,G2i,G2j,G3i,G3j]';

vP=jacobian(P,q)*vq+jacobian(P,G3i)*vG3i+jacobian(P,G3j)*vG3j;
aP=jacobian(vP,q)*vq+jacobian(vP,vq)*aq+jacobian(vP,G3i)*vG3i+jacobian(vP,G3j)*vG3j+jacobian(vP,vG3i)*aG3i+jacobian(vP,vG3j)*aG3j;

vG=jacobian(G,q)*vq+jacobian(G,G3i)*vG3i+jacobian(G,G3j)*vG3j;
aG=jacobian(vG,q)*vq+jacobian(vG,vq)*aq+jacobian(vG,G3i)*vG3i+jacobian(vG,G3j)*vG3j+jacobian(vG,vG3i)*aG3i+jacobian(vG,vG3j)*aG3j;

%     % Kinematics (body based):
%     G3i=G_list(5); G3j=G_list(6); vG3i=vG_list(5); vG3j=vG_list(6); aG3i=aG_list(5); aG3j=aG_list(6);
%     q1=q_list(1); q2=q_list(2); q3=q_list(3);
%     vq1=vq_list(1); vq2=vq_list(2); vq3=vq_list(3);
%     aq1=aq_list(1); aq2=aq_list(2); aq3=aq_list(3);
%     
%     P3i=G3i+(L3-d3)*sin(q3); P3j=G3j+(L3-d3)*cos(q3);
%     P2i=G3i-d3*sin(q3); P2j=G3j-d3*cos(q3);
%     P1i=P2i-L2*sin(q2); P1j=P2j-L2*cos(q2);
%     P0i=P1i-L1*sin(q1); P0j=P1j-L1*cos(q1);
%     P=[P0i,P0j,P1i,P1j,P2i,P2j,P3i,P3j]';
% 
%     G2i=P2i-(L2-d2)*sin(q2); G2j=P2j-(L2-d2)*cos(q2);
%     G1i=P1i-(L1-d1)*sin(q1); G1j=P1j-(L1-d1)*cos(q1);
%     G=[G1i,G1j,G2i,G2j,G3i,G3j]';
%     
%     vP0i=vG3i - vq1*cos(conj(q1))*conj(L1) - vq2*cos(conj(q2))*conj(L2) - vq3*cos(conj(q3))*conj(d3);
%     vP0j=vG3j + vq1*sin(conj(q1))*conj(L1) + vq2*sin(conj(q2))*conj(L2) + vq3*sin(conj(q3))*conj(d3);
%     vP0=[vP0i,vP0j]';
% 
%     aP0i=sin(conj(q1))*conj(L1)*vq1^2 + sin(conj(q2))*conj(L2)*vq2^2 + sin(conj(q3))*conj(d3)*vq3^2 + aG3i - aq1*cos(conj(q1))*conj(L1) - aq2*cos(conj(q2))*conj(L2) - aq3*cos(conj(q3))*conj(d3);
%     aP0j=cos(conj(q1))*conj(L1)*vq1^2 + cos(conj(q2))*conj(L2)*vq2^2 + cos(conj(q3))*conj(d3)*vq3^2 + aG3j + aq1*sin(conj(q1))*conj(L1) + aq2*sin(conj(q2))*conj(L2) + aq3*sin(conj(q3))*conj(d3);
%     aP0=[aP0i,aP0j]';
%     
%     vG1i=vG3i - vq1*cos(conj(q1))*(conj(L1) - conj(d1)) - vq2*cos(conj(q2))*conj(L2) - vq3*cos(conj(q3))*conj(d3);
%     vG1j=vG3j + vq1*sin(conj(q1))*(conj(L1) - conj(d1)) + vq2*sin(conj(q2))*conj(L2) + vq3*sin(conj(q3))*conj(d3);
%     vG2i=vG3i - vq2*cos(conj(q2))*(conj(L2) - conj(d2)) - vq3*cos(conj(q3))*conj(d3);
%     vG2j=vG3j + vq2*sin(conj(q2))*(conj(L2) - conj(d2)) + vq3*sin(conj(q3))*conj(d3);
%     vG=[vG1i,vG1j,vG2i,vG2j,vG3i,vG3j]';
%     
%     aG1i=sin(conj(q1))*(conj(L1) - conj(d1))*vq1^2 + sin(conj(q2))*conj(L2)*vq2^2 + sin(conj(q3))*conj(d3)*vq3^2 + aG3i - aq1*cos(conj(q1))*(conj(L1) - conj(d1)) - aq2*cos(conj(q2))*conj(L2) - aq3*cos(conj(q3))*conj(d3);
%     aG1j=cos(conj(q1))*(conj(L1) - conj(d1))*vq1^2 + cos(conj(q2))*conj(L2)*vq2^2 + cos(conj(q3))*conj(d3)*vq3^2 + aG3j + aq1*sin(conj(q1))*(conj(L1) - conj(d1)) + aq2*sin(conj(q2))*conj(L2) + aq3*sin(conj(q3))*conj(d3);
%     aG2i=sin(conj(q2))*(conj(L2) - conj(d2))*vq2^2 + sin(conj(q3))*conj(d3)*vq3^2 + aG3i - aq2*cos(conj(q2))*(conj(L2) - conj(d2)) - aq3*cos(conj(q3))*conj(d3);
%     aG2j=cos(conj(q2))*(conj(L2) - conj(d2))*vq2^2 + cos(conj(q3))*conj(d3)*vq3^2 + aG3j + aq2*sin(conj(q2))*(conj(L2) - conj(d2)) + aq3*sin(conj(q3))*conj(d3);
%     aG=[aG1i,aG1j,aG2i,aG2j,aG3i,aG3j]';
%     
%     ceq=[ceq,(P_list-P)',(G_list-G)',(vP0_list-vP0)',(aP0_list-aP0)',(vG_list-vG)',(aG_list-aG)'];


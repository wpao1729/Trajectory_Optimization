function [c,ceq] = step_strike_nlcon (var_list,n,m_list,I_list,g,L_list,d_list)
% Nonlinear conditions.
% P,G,vG,aG are symbolic expressions dependent on q,vq,aq.


ceq=[]; c=0;

q_mat=reshape(var_list(1:(5*n)),5,n);
vq_mat=reshape(var_list((5*n+1):(10*n)),5,n);
aq_mat=reshape(var_list((10*n+1):(15*n)),5,n);
P_mat=reshape(var_list((15*n+1):(25*n)),10,n);
G_mat=reshape(var_list((25*n+1):(35*n)),10,n);
vG_mat=reshape(var_list((35*n+1):(45*n)),10,n);
aG_mat=reshape(var_list((45*n+1):(55*n)),10,n);
u_mat=reshape(var_list((55*n+1):(59*n)),4,n);

L1=L_list(1);L2=L_list(2);L3=L_list(3);L4=L_list(4);L5=L_list(5);
d1=d_list(1);d2=d_list(2);d3=d_list(3);d4=d_list(4);d5=d_list(5);

% single stance kinematics and dynamics
for i=1:(n-1) % at every knot point
    
    q_list=q_mat(:,i); vq_list=vq_mat(:,i); aq_list=aq_mat(:,i); P_list=P_mat(:,i); G_list=G_mat(:,i); vG_list=vG_mat(:,i); aG_list=aG_mat(:,i); u_list=[0;u_mat(:,i)];
    
    % Dynamics:
    P0i=0; P0j=0; Pu_list=[P0i;P0j;P_list(1:4);P_list(3:4);P_list(7:8)];
    Pui_list=Pu_list(1:2:10); Puj_list=Pu_list(2:2:10); Gi_list=G_list(1:2:10); Gj_list=G_list(2:2:10); aGi_list=aG_list(1:2:10); aGj_list=aG_list(2:2:10); 
    for k=1:length(u_list)
        Pui=Pui_list(k); Puj=Puj_list(k);
        torqueG_list=-m_list*g.*(Gi_list-Pui);
        accel_list=m_list.*((Gi_list-Pui).*aGj_list-(Gj_list-Puj).*aGi_list)-I_list.*aq_list;
        eq_dynamics=u_list(k)+sum(torqueG_list(k:5))-sum(accel_list(k:5));
        ceq=[ceq,eq_dynamics];
    end

    % Kinematics:
    % !!! PUT KINEMATICS EXPLICITLY HERE!!!
    q1=q_list(1); q2=q_list(2); q3=q_list(3); q4=q_list(4); q5=q_list(5);
    vq1=vq_list(1); vq2=vq_list(2); vq3=vq_list(3); vq4=vq_list(4); vq5=vq_list(5);
    aq1=aq_list(1); aq2=aq_list(2); aq3=aq_list(3); aq4=aq_list(4); aq5=aq_list(5);
    
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

    vG1i=(vq1*cos((q1)))*d1;
    vG1j=-(vq1*sin((q1)))*d1;
    vG2i=(vq1*cos((q1)))*L1 + (vq2*cos((q2)))*d2;
    vG2j=- (vq1*sin((q1)))*L1 - (vq2*sin((q2)))*d2;
    vG3i=(vq1*cos((q1)))*L1 + (vq2*cos((q2)))*L2 + (vq3*cos((q3)))*d3;
    vG3j=- (vq1*sin((q1)))*L1 - (vq2*sin((q2)))*L2 - (vq3*sin((q3)))*d3;
    vG4i=(vq1*cos((q1)))*L1 + (vq2*cos((q2)))*L2 - (vq4*cos((q4)))*d4;
    vG4j=(vq4*sin((q4)))*L1 - (vq2*sin((q2)))*L2 - (vq1*sin((q1)))*d4;
    vG5i=(vq1*cos((q1)))*L1 + (vq2*cos((q2)))*L2 - (vq4*cos((q4)))*L4 - (vq5*cos((q5)))*d5;
    vG5j=(vq4*sin((q4)))*L1 - (vq2*sin((q2)))*L2 - (vq1*sin((q1)))*L4 + (vq5*sin((q5)))*d5;
    vG=[vG1i,vG1j,vG2i,vG2j,vG3i,vG3j,vG4i,vG4j,vG5i,vG5j]';
    
    aG1i=- sin((q1))*(d1)*vq1^2 + aq1*cos((q1))*(d1);
    aG1j=- cos((q1))*(d1)*vq1^2 - aq1*sin((q1))*(d1);
    aG2i=- sin((q1))*(L1)*vq1^2 - sin((q2))*(d2)*vq2^2 + aq1*cos((q1))*(L1) + aq2*cos((q2))*(d2);
    aG2j=- cos((q1))*(L1)*vq1^2 - cos((q2))*(d2)*vq2^2 - aq1*sin((q1))*(L1) - aq2*sin((q2))*(d2);
    aG3i=- sin((q1))*(L1)*vq1^2 - sin((q2))*(L2)*vq2^2 - sin((q3))*(d3)*vq3^2 + aq1*cos((q1))*(L1) + aq2*cos((q2))*(L2) + aq3*cos((q3))*(d3);
    aG3j=- cos((q1))*(L1)*vq1^2 - cos((q2))*(L2)*vq2^2 - cos((q3))*(d3)*vq3^2 - aq1*sin((q1))*(L1) - aq2*sin((q2))*(L2) - aq3*sin((q3))*(d3);
    aG4i=- sin((q1))*(L1)*vq1^2 - sin((q2))*(L2)*vq2^2 + sin((q4))*(d4)*vq4^2 + aq1*cos((q1))*(L1) + aq2*cos((q2))*(L2) - aq4*cos((q4))*(d4);
    aG4j=- cos((q1))*(L1)*vq1^2 - cos((q2))*(L2)*vq2^2 + cos((q4))*(d4)*vq4^2 - aq1*sin((q1))*(L1) - aq2*sin((q2))*(L2) + aq4*sin((q4))*(d4);
    aG5i=- sin((q1))*(L1)*vq1^2 - sin((q2))*(L2)*vq2^2 + sin((q4))*(L4)*vq4^2 + sin((q5))*(d5)*vq5^2 + aq1*cos((q1))*(L1) + aq2*cos((q2))*(L2) - aq4*cos((q4))*(L4) - aq5*cos((q5))*(d5);
    aG5j=- cos((q1))*(L1)*vq1^2 - cos((q2))*(L2)*vq2^2 + cos((q4))*(L4)*vq4^2 + cos((q5))*(d5)*vq5^2 - aq1*sin((q1))*(L1) - aq2*sin((q2))*(L2) + aq4*sin((q4))*(L4) + aq5*sin((q5))*(d5);
    aG=[aG1i,aG1j,aG2i,aG2j,aG3i,aG3j,aG4i,aG4j,aG5i,aG5j]';

    ceq=[ceq,(P_list-P)',(G_list-G)',(vG_list-vG)',(aG_list-aG)'];
    
end


%strike dynamics (conservation of angular momentum)
%conservation of angular momentum at front stance
% Pu_list1=P_mat([1:4,3:4,7:10],n-1); Pui_list1=Pu_list1(1:2:10); Puj_list1=Pu_list1(2:2:10);
% Pu_list2=P_mat([1:4,3:4,7:10],n); Pui_list2=Pu_list2(1:2:10); Puj_list2=Pu_list2(2:2:10);
% for k=5 %1:length(u_list)
%     Pui1=Pui_list1(k); Puj1=Puj_list1(k); Pui2=Pui_list2(k); Puj2=Puj_list2(k);
%     momentum_list1=m_list.*((G_mat(1:2:9,(n-1))-Pui1).*vG_mat(2:2:10,(n-1))-(G_mat(2:2:10,(n-1))-Puj1).*vG_mat(1:2:9,(n-1)))-I_list.*vq_mat(:,(n-1));
%     momentum_list2=m_list.*((G_mat(1:2:9,n)-Pui2).*vG_mat(2:2:10,n)-(G_mat(2:2:10,n)-Puj2).*vG_mat(1:2:9,n))-I_list.*vq_mat(:,n);
%     eq_strike=sum(momentum_list1(1:k))-sum(momentum_list2(1:k));
%     ceq=[ceq,eq_strike];
% end
% conservation of angular momentum at back stance
% Pui1=0; Puj1=0; Pui2=0; Puj2=0;
% momentum_list1=m_list.*((G_mat(1:2:9,(n-1))-Pui1).*vG_mat(2:2:10,(n-1))-(G_mat(2:2:10,(n-1))-Puj1).*vG_mat(1:2:9,(n-1)))-I_list.*vq_mat(:,(n-1));
% momentum_list2=m_list.*((G_mat(1:2:9,n)-Pui2).*vG_mat(2:2:10,n)-(G_mat(2:2:10,n)-Puj2).*vG_mat(1:2:9,n))-I_list.*vq_mat(:,n);
% eq_strike=sum(momentum_list1(1:5))-sum(momentum_list2(1:5));
% ceq=[ceq,eq_strike];
% conservation of momentum in i direction
% eq_strike=sum(m_list.*G_mat(1:2:9,n-1))-sum(m_list.*G_mat(1:2:9,n));
% ceq=[ceq,eq_strike];

% Second step dynamics and kinematics (single stance)
% Dynamics:
P_list=P_mat(:,n);
Pu_list=[P_list(1:4);P_list(3:4);P_list(7:10)];
Pui_list=Pu_list(1:2:10); Puj_list=Pu_list(2:2:10);
u_list=[u_mat(:,n);0];
for k=1:length(u_list)
    Pui=Pui_list(k); Puj=Puj_list(k);
    torqueG_list=-m_list*g.*(G_mat(1:2:9,n)-Pui);
    accel_list=m_list.*((G_mat(1:2:9,n)-Pui).*aG_mat(2:2:10,n)-(G_mat(2:2:10)-Puj).*aG_mat(1:2:9,n))-I_list.*aq_mat(:,n);
    eq_dynamics=-u_list(k)+sum(torqueG_list(1:k))-sum(accel_list(1:k));
    ceq=[ceq,eq_dynamics];
end

% Kinematics:
q_list=q_mat(:,n); vq_list=vq_mat(:,n); aq_list=aq_mat(:,n); P_list=P_mat(:,n); G_list=G_mat(:,n); vG_list=vG_mat(:,n); aG_list=aG_mat(:,n);
q1=q_list(1); q2=q_list(2); q3=q_list(3); q4=q_list(4); q5=q_list(5);
vq1=vq_list(1); vq2=vq_list(2); vq3=vq_list(3); vq4=vq_list(4); vq5=vq_list(5);
aq1=aq_list(1); aq2=aq_list(2); aq3=aq_list(3); aq4=aq_list(4); aq5=aq_list(5);

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

vG1i=vq4*cos((q4))*(L4) - vq2*cos((q2))*(L2) - vq1*cos((q1))*((L1) - (d1)) + vq5*cos((q5))*(L5);
vG1j=vq1*sin((q1))*((L1) - (d1)) + vq2*sin((q2))*(L2) - vq4*sin((q4))*(L4) - vq5*sin((q5))*(L5);
vG2i=vq4*cos((q4))*(L4) - vq2*cos((q2))*((L2) - (d2)) + vq5*cos((q5))*(L5);
vG2j=vq2*sin((q2))*((L2) - (d2)) - vq4*sin((q4))*(L4) - vq5*sin((q5))*(L5);
vG3i=vq4*cos((q4))*(L4) + vq5*cos((q5))*(L5) + vq3*cos((q3))*(d3);
vG3j=- vq4*sin((q4))*(L4) - vq5*sin((q5))*(L5) - vq3*sin((q3))*(d3);
vG4i=vq4*cos((q4))*((L4) - (d4)) + vq5*cos((q5))*(L5);
vG4j=- vq4*sin((q4))*((L4) - (d4)) - vq5*sin((q5))*(L5);
vG5i=vq5*cos((q5))*((L5) - (d5));
vG5j=-vq5*sin((q5))*((L5) - (d5));
vG=[vG1i,vG1j,vG2i,vG2j,vG3i,vG3j,vG4i,vG4j,vG5i,vG5j]';

aG1i=sin((q1))*((L1) - (d1))*vq1^2 + sin((q2))*(L2)*vq2^2 - sin((q4))*(L4)*vq4^2 - sin((q5))*(L5)*vq5^2 - aq1*cos((q1))*((L1) - (d1)) - aq2*cos((q2))*(L2) + aq4*cos((q4))*(L4) + aq5*cos((q5))*(L5);
aG1j=cos((q1))*((L1) - (d1))*vq1^2 + cos((q2))*(L2)*vq2^2 - cos((q4))*(L4)*vq4^2 - cos((q5))*(L5)*vq5^2 + aq1*sin((q1))*((L1) - (d1)) + aq2*sin((q2))*(L2) - aq4*sin((q4))*(L4) - aq5*sin((q5))*(L5);
aG2i=sin((q2))*((L2) - (d2))*vq2^2 - sin((q4))*(L4)*vq4^2 - sin((q5))*(L5)*vq5^2 - aq2*cos((q2))*((L2) - (d2)) + aq4*cos((q4))*(L4) + aq5*cos((q5))*(L5);
aG2j=cos((q2))*((L2) - (d2))*vq2^2 - cos((q4))*(L4)*vq4^2 - cos((q5))*(L5)*vq5^2 + aq2*sin((q2))*((L2) - (d2)) - aq4*sin((q4))*(L4) - aq5*sin((q5))*(L5);
aG3i=- sin((q3))*(d3)*vq3^2 - sin((q4))*(L4)*vq4^2 - sin((q5))*(L5)*vq5^2 + aq4*cos((q4))*(L4) + aq5*cos((q5))*(L5) + aq3*cos((q3))*(d3);
aG3j=- cos((q3))*(d3)*vq3^2 - cos((q4))*(L4)*vq4^2 - cos((q5))*(L5)*vq5^2 - aq4*sin((q4))*(L4) - aq5*sin((q5))*(L5) - aq3*sin((q3))*(d3);
aG4i=- sin((q4))*((L4) - (d4))*vq4^2 - sin((q5))*(L5)*vq5^2 + aq4*cos((q4))*((L4) - (d4)) + aq5*cos((q5))*(L5);
aG4j=- cos((q4))*((L4) - (d4))*vq4^2 - cos((q5))*(L5)*vq5^2 - aq4*sin((q4))*((L4) - (d4)) - aq5*sin((q5))*(L5);
aG5i=- sin((q5))*((L5) - (d5))*vq5^2 + aq5*cos((q5))*((L5) - (d5));
aG5j=- cos((q5))*((L5) - (d5))*vq5^2 - aq5*sin((q5))*((L5) - (d5));

aG=[aG1i,aG1j,aG2i,aG2j,aG3i,aG3j,aG4i,aG4j,aG5i,aG5j]';

ceq=[ceq,(P_list-P)',(G_list-G)',(vG_list-vG)',(aG_list-aG)'];


end
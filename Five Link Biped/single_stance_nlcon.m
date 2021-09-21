% function [c,ceq] = single_stance_nlcon (var_list,n,m_list,I_list,g,P,G,vG,aG)
% % Nonlinear conditions.
% % P,G,vG,aG are symbolic expressions dependent on q,vq,aq.
% 
% 
% ceq=[]; c=0;
% 
% q_mat=reshape(var_list(1:(5*n)),5,n);
% vq_mat=reshape(var_list((5*n+1):(10*n)),5,n);
% aq_mat=reshape(var_list((10*n+1):(15*n)),5,n);
% P_mat=reshape(var_list((15*n+1):(25*n)),10,n);
% G_mat=reshape(var_list((25*n+1):(35*n)),10,n);
% vG_mat=reshape(var_list((35*n+1):(45*n)),10,n);
% aG_mat=reshape(var_list((45*n+1):(55*n)),10,n);
% u_mat=reshape(var_list((55*n+1):(60*n)),5,n);
% 
% for i=1:n % at every knot point
%     
%     q_list=q_mat(:,i); vq_list=vq_mat(:,i); aq_list=aq_mat(:,i); P_list=P_mat(:,i); G_list=G_mat(:,i); vG_list=vG_mat(:,i); aG_list=aG_mat(:,i); u_list=u_mat(:,i);
%     
%     % Dynamics:
%     P0i=0; P0j=0; P2_list=[P0i;P0j;P_list(1:4);P_list(3:4);P_list(7:8)];
%     P2i_list=P2_list(1:2:10); P2j_list=P2_list(2:2:10); Gi_list=G_list(1:2:10); Gj_list=G_list(2:2:10); aGi_list=aG_list(1:2:10); aGj_list=aG_list(2:2:10); 
%     for k=1:length(u_list)
%         P2i=P2i_list(k); P2j=P2j_list(k);
%         torqueG_list=-m_list*g.*(Gi_list-P2i);
%         accel_list=m_list.*((Gi_list-P2i).*aGj_list-(Gj_list-P2j).*aGi_list)-I_list.*aq_list;
%         eq_dynamics=u_list(k)+sum(torqueG_list(k:5))-sum(accel_list(k:5));
%         ceq=[ceq,eq_dynamics];
%     end
% 
%     % Kinematics:
%     % !!! PUT KINEMATICS EXPLICITLY HERE!!!
%     
%     P_knm_list=P; G_knm_list=G; vG_knm_list=vG; aG_knm_list=aG; %_knm_list contains values calculated by kinematics
%     for k=1:length(q_list)
%         % Replace symbolic variables in the equations with corresponding values
%         P_knm_list=subs(P_knm_list,['q',num2str(k)],q_list(k));
%         G_knm_list=subs(G_knm_list,['q',num2str(k)],q_list(k));
%         vG_knm_list=subs(vG_knm_list,['q',num2str(k)],q_list(k));
%         vG_knm_list=subs(vG_knm_list,['vq',num2str(k)],vq_list(k));
%         aG_knm_list=subs(aG_knm_list,['q',num2str(k)],q_list(k));
%         aG_knm_list=subs(aG_knm_list,['vq',num2str(k)],vq_list(k));
%         aG_knm_list=subs(aG_knm_list,['aq',num2str(k)],aq_list(k));
%     end
%     P_knm_list=double(P_knm_list); G_knm_list=double(G_knm_list); vG_knm_list=double(vG_knm_list); aG_knm_list=double(aG_knm_list);
%     ceq=[ceq,(P_list-P_knm_list)',(G_list-G_knm_list)',(vG_list-vG_knm_list)',(aG_list-aG_knm_list)'];
%     
% end
% 
% 
% end

function [c,ceq] = single_stance_nlcon (var_list,n,m_list,I_list,g,L_list,d_list)
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
u_mat=reshape(var_list((55*n+1):(60*n)),5,n);

L1=L_list(1);L2=L_list(2);L3=L_list(3);L4=L_list(4);L5=L_list(5);
d1=d_list(1);d2=d_list(2);d3=d_list(3);d4=d_list(4);d5=d_list(5);

for i=1:n % at every knot point
    
    q_list=q_mat(:,i); vq_list=vq_mat(:,i); aq_list=aq_mat(:,i); P_list=P_mat(:,i); G_list=G_mat(:,i); vG_list=vG_mat(:,i); aG_list=aG_mat(:,i); u_list=u_mat(:,i);
    
    % Dynamics:
    P0i=0; P0j=0; P2_list=[P0i;P0j;P_list(1:4);P_list(3:4);P_list(7:8)];
    P2i_list=P2_list(1:2:10); P2j_list=P2_list(2:2:10); Gi_list=G_list(1:2:10); Gj_list=G_list(2:2:10); aGi_list=aG_list(1:2:10); aGj_list=aG_list(2:2:10); 
    for k=1:length(u_list)
        P2i=P2i_list(k); P2j=P2j_list(k);
        torqueG_list=-m_list*g.*(Gi_list-P2i);
        accel_list=m_list.*((Gi_list-P2i).*aGj_list-(Gj_list-P2j).*aGi_list)-I_list.*aq_list;
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


end
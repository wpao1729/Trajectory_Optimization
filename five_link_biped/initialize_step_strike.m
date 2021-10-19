function var_list_guess=initialize_step_strike(n,t_list,qi_list,qf_list,L_list,d_list,m_list,I_list,g)
    
L1=L_list(1);L2=L_list(2);L3=L_list(3);L4=L_list(4);L5=L_list(5);
d1=d_list(1);d2=d_list(2);d3=d_list(3);d4=d_list(4);d5=d_list(5);

var_list_guess=zeros(1,59*n);
q_mat=zeros(5,n); vq_mat=zeros(5,n); aq_mat=zeros(5,n);
tf=t_list(end);
for k=1:5
    a=-6*(qf_list(k)-qi_list(k))/(tf^3); 
    q_mat(k,:)=a/3*(t_list.^3)-a*tf/2*(t_list.^2)+qi_list(k);
    vq_mat(k,:)=a*(t_list.^2)-a*tf*t_list;
    aq_mat(k,:)=2*a*t_list-a*tf;
end
var_list_guess(1:(5*n))=reshape(q_mat,1,5*n);
var_list_guess((5*n+1):(10*n))=reshape(vq_mat,1,5*n);
var_list_guess((10*n+1):(15*n))=reshape(aq_mat,1,5*n);

P_mat=zeros(10,n); G_mat=zeros(10,n); vG_mat=zeros(10,n); aG_mat=zeros(10,n); u_mat=zeros(4,n);
for i=1:n % at every knot point

    q_list=q_mat(:,i); vq_list=vq_mat(:,i); aq_list=aq_mat(:,i);

    % Kinematics:
    q1=q_list(1); q2=q_list(2); q3=q_list(3); q4=q_list(4); q5=q_list(5);
    vq1=vq_list(1); vq2=vq_list(2); vq3=vq_list(3); vq4=vq_list(4); vq5=vq_list(5);
    aq1=aq_list(1); aq2=aq_list(2); aq3=aq_list(3); aq4=aq_list(4); aq5=aq_list(5);

    P1i=L1*sin(q1); P1j=L1*cos(q1);
    P2i=P1i+L2*sin(q2); P2j=P1j+L2*cos(q2);
    P3i=P2i+L3*sin(q3); P3j=P2j+L3*cos(q3);
    P4i=P2i-L4*sin(q4); P4j=P2j-L4*cos(q4);
    P5i=P4i-L5*sin(q5); P5j=P4j-L5*cos(q5);
    P_mat(:,i)=[P1i,P1j,P2i,P2j,P3i,P3j,P4i,P4j,P5i,P5j]';

    G1i=d1*sin(q1); G1j=d1*cos(q1);
    G2i=P1i+d2*sin(q2); G2j=P1j+d2*cos(q2);
    G3i=P2i+d3*sin(q3); G3j=P2j+d3*cos(q3);
    G4i=P2i-d4*sin(q4); G4j=P2j-d4*cos(q4);
    G5i=P4i-d5*sin(q5); G5j=P4j-d5*cos(q5);
    G_mat(:,i)=[G1i,G1j,G2i,G2j,G3i,G3j,G4i,G4j,G5i,G5j]';

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
    vG_mat(:,i)=[vG1i,vG1j,vG2i,vG2j,vG3i,vG3j,vG4i,vG4j,vG5i,vG5j]';

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
    aG_mat(:,i)=[aG1i,aG1j,aG2i,aG2j,aG3i,aG3j,aG4i,aG4j,aG5i,aG5j]';
    
    % Dynamics:
    P0i=0; P0j=0; P2_list=[P0i;P0j;P_mat(1:4,i);P_mat(3:4,i);P_mat(7:8,i)];
    P2i_list=P2_list(1:2:10); P2j_list=P2_list(2:2:10); Gi_list=G_mat(1:2:10,i); Gj_list=G_mat(2:2:10,i); aGi_list=aG_mat(1:2:10,i); aGj_list=aG_mat(2:2:10,i); 
    u_list=zeros(5,1);
    for k=1:5
        P2i=P2i_list(k); P2j=P2j_list(k);
        torqueG_list=-m_list*g.*(Gi_list-P2i);
        accel_list=m_list.*((Gi_list-P2i).*aGj_list-(Gj_list-P2j).*aGi_list)-I_list.*aq_list;
        u_list(k)=-sum(torqueG_list(k:5))+sum(accel_list(k:5));
    end
    u_mat(:,i)=u_list(2:5);
end

var_list_guess((15*n+1):(25*n))=reshape(P_mat,1,10*n);
var_list_guess((25*n+1):(35*n))=reshape(G_mat,1,10*n);
var_list_guess((35*n+1):(45*n))=reshape(vG_mat,1,10*n);
var_list_guess((45*n+1):(55*n))=reshape(aG_mat,1,10*n);
var_list_guess((55*n+1):(59*n))=reshape(u_mat,1,4*n);

end
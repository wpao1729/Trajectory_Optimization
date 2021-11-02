function var_list_guess=initialize_running(n,qi_list,qf_list,L_list,d_list,m_list,I_list,g,d_run,h_run,t_list,iswalking)

var_list_guess=zeros(1,75*n+2);

L1=L_list(1);L2=L_list(2);L3=L_list(3);L4=L_list(4);L5=L_list(5);
d1=d_list(1);d2=d_list(2);d3=d_list(3);d4=d_list(4);d5=d_list(5);

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

% Initialize (av)G3ij
G3i_i=L1*sin(q_mat(1,1))+L2*sin(q_mat(2,1))+d3*sin(q_mat(3,1));
G3j_i=L1*cos(q_mat(1,1))+L2*cos(q_mat(2,1))+d3*cos(q_mat(3,1));
G3i_f=d_run+L5*sin(q_mat(5,end))+L4*sin(q_mat(4,end))+d3*sin(q_mat(3,end));
G3j_f=h_run+L5*cos(q_mat(5,end))+L4*cos(q_mat(4,end))+d3*cos(q_mat(3,end));
a_G3i=-6*(G3i_f-G3i_i)/(tf^3); 
G3i_list=a_G3i/3*(t_list.^3)-a_G3i*tf/2*(t_list.^2)+G3i_i;
vG3i_list=a_G3i*(t_list.^2)-a_G3i*tf*t_list;
aG3i_list=2*a_G3i*t_list-a_G3i*tf;
a_G3j=-6*(G3j_f-G3j_i)/(tf^3); 
G3j_list=a_G3j/3*(t_list.^3)-a_G3j*tf/2*(t_list.^2)+G3j_i;
vG3j_list=a_G3j*(t_list.^2)-a_G3j*tf*t_list;
aG3j_list=2*a_G3j*t_list-a_G3j*tf;


P_mat=zeros(12,n); vP0_mat=zeros(2,n); aP0_mat=zeros(2,n); vP5_mat=zeros(2,n); aP5_mat=zeros(2,n); G_mat=zeros(10,n); vG_mat=zeros(10,n); aG_mat=zeros(10,n); GRF_mat=zeros(4,n); u_mat=zeros(6,n);
for i=1:n % at every knot point

    q_list=q_mat(:,i); vq_list=vq_mat(:,i); aq_list=aq_mat(:,i);

    % Kinematics:
    G3i=G3i_list(i); G3j=G3j_list(i); vG3i=vG3i_list(i); vG3j=vG3j_list(i); aG3i=aG3i_list(i); aG3j=aG3j_list(i);
    
    q1=q_list(1); q2=q_list(2); q3=q_list(3); q4=q_list(4); q5=q_list(5);
    vq1=vq_list(1); vq2=vq_list(2); vq3=vq_list(3); vq4=vq_list(4); vq5=vq_list(5);
    aq1=aq_list(1); aq2=aq_list(2); aq3=aq_list(3); aq4=aq_list(4); aq5=aq_list(5);
    
    P3i=G3i+(L3-d3)*sin(q3); P3j=G3j+(L3-d3)*cos(q3);
    P2i=G3i-d3*sin(q3); P2j=G3j-d3*cos(q3);
    P1i=P2i-L2*sin(q2); P1j=P2j-L2*cos(q2);
    P0i=P1i-L1*sin(q1); P0j=P1j-L1*cos(q1);
    P4i=P2i-L4*sin(q4); P4j=P2j-L4*cos(q4);
    P5i=P4i-L5*sin(q5); P5j=P4j-L5*cos(q5);
    P_mat(:,i)=[P0i,P0j,P1i,P1j,P2i,P2j,P3i,P3j,P4i,P4j,P5i,P5j]';

    G2i=P2i-d2*sin(q2); G2j=P2j-d2*cos(q2);
    G1i=P1i-d1*sin(q1); G1j=P1j-d1*cos(q1);
    G4i=P2i-d4*sin(q4); G4j=P2j-d4*cos(q4);
    G5i=P4i-d5*sin(q5); G5j=P4j-d5*cos(q5);
    G_mat(:,i)=[G1i,G1j,G2i,G2j,G3i,G3j,G4i,G4j,G5i,G5j]';
    
    vP0i=vG3i - vq1*cos(conj(q1))*conj(L1) - vq2*cos(conj(q2))*conj(L2) - vq3*cos(conj(q3))*conj(d3);
    vP0j=vG3j + vq1*sin(conj(q1))*conj(L1) + vq2*sin(conj(q2))*conj(L2) + vq3*sin(conj(q3))*conj(d3);
    vP0_mat(:,i)=[vP0i,vP0j]';
    vP5i=vG3i - vq4*cos(conj(q4))*conj(L4) - vq5*cos(conj(q5))*conj(L5) - vq3*cos(conj(q3))*conj(d3);
    vP5j=vG3j + vq4*sin(conj(q4))*conj(L4) + vq5*sin(conj(q5))*conj(L5) + vq3*sin(conj(q3))*conj(d3);
    vP5_mat(:,i)=[vP5i,vP5j]';
    
    aP0i=sin(conj(q1))*conj(L1)*vq1^2 + sin(conj(q2))*conj(L2)*vq2^2 + sin(conj(q3))*conj(d3)*vq3^2 + aG3i - aq1*cos(conj(q1))*conj(L1) - aq2*cos(conj(q2))*conj(L2) - aq3*cos(conj(q3))*conj(d3);
    aP0j=cos(conj(q1))*conj(L1)*vq1^2 + cos(conj(q2))*conj(L2)*vq2^2 + cos(conj(q3))*conj(d3)*vq3^2 + aG3j + aq1*sin(conj(q1))*conj(L1) + aq2*sin(conj(q2))*conj(L2) + aq3*sin(conj(q3))*conj(d3);
    aP0_mat(:,i)=[aP0i,aP0j]';
    aP5i=sin(conj(q3))*conj(d3)*vq3^2 + sin(conj(q4))*conj(L4)*vq4^2 + sin(conj(q5))*conj(L5)*vq5^2 + aG3i - aq4*cos(conj(q4))*conj(L4) - aq5*cos(conj(q5))*conj(L5) - aq3*cos(conj(q3))*conj(d3);
    aP5j=cos(conj(q3))*conj(d3)*vq3^2 + cos(conj(q4))*conj(L4)*vq4^2 + cos(conj(q5))*conj(L5)*vq5^2 + aG3j + aq4*sin(conj(q4))*conj(L4) + aq5*sin(conj(q5))*conj(L5) + aq3*sin(conj(q3))*conj(d3);
    aP5_mat(:,i)=[aP5i,aP5j]';
    
    vG1i=vG3i - vq2*cos(conj(q2))*conj(L2) - vq1*cos(conj(q1))*conj(d1) - vq3*cos(conj(q3))*conj(d3);
    vG1j=vG3j + vq2*sin(conj(q2))*conj(L2) + vq1*sin(conj(q1))*conj(d1) + vq3*sin(conj(q3))*conj(d3);
    vG2i=vG3i - vq2*cos(conj(q2))*conj(d2) - vq3*cos(conj(q3))*conj(d3);
    vG2j=vG3j + vq2*sin(conj(q2))*conj(d2) + vq3*sin(conj(q3))*conj(d3);
    vG4i=vG3i - vq3*cos(conj(q3))*conj(d3) - vq4*cos(conj(q4))*conj(d4);
    vG4j=vG3j + vq3*sin(conj(q3))*conj(d3) + vq4*sin(conj(q4))*conj(d4);
    vG5i=vG3i - vq4*cos(conj(q4))*conj(L4) - vq3*cos(conj(q3))*conj(d3) - vq5*cos(conj(q5))*conj(d5);
    vG5j=vG3j + vq4*sin(conj(q4))*conj(L4) + vq3*sin(conj(q3))*conj(d3) + vq5*sin(conj(q5))*conj(d5);
    vG_mat(:,i)=[vG1i,vG1j,vG2i,vG2j,vG3i,vG3j,vG4i,vG4j,vG5i,vG5j]';
    
    aG1i=sin(conj(q1))*conj(d1)*vq1^2 + sin(conj(q2))*conj(L2)*vq2^2 + sin(conj(q3))*conj(d3)*vq3^2 + aG3i - aq2*cos(conj(q2))*conj(L2) - aq1*cos(conj(q1))*conj(d1) - aq3*cos(conj(q3))*conj(d3);
    aG1j=cos(conj(q1))*conj(d1)*vq1^2 + cos(conj(q2))*conj(L2)*vq2^2 + cos(conj(q3))*conj(d3)*vq3^2 + aG3j + aq2*sin(conj(q2))*conj(L2) + aq1*sin(conj(q1))*conj(d1) + aq3*sin(conj(q3))*conj(d3);
    aG2i=sin(conj(q2))*conj(d2)*vq2^2 + sin(conj(q3))*conj(d3)*vq3^2 + aG3i - aq2*cos(conj(q2))*conj(d2) - aq3*cos(conj(q3))*conj(d3);
    aG2j=cos(conj(q2))*conj(d2)*vq2^2 + cos(conj(q3))*conj(d3)*vq3^2 + aG3j + aq2*sin(conj(q2))*conj(d2) + aq3*sin(conj(q3))*conj(d3);
    aG4i=sin(conj(q3))*conj(d3)*vq3^2 + sin(conj(q4))*conj(d4)*vq4^2 + aG3i - aq3*cos(conj(q3))*conj(d3) - aq4*cos(conj(q4))*conj(d4);
    aG4j=cos(conj(q3))*conj(d3)*vq3^2 + cos(conj(q4))*conj(d4)*vq4^2 + aG3j + aq3*sin(conj(q3))*conj(d3) + aq4*sin(conj(q4))*conj(d4);
    aG5i=sin(conj(q3))*conj(d3)*vq3^2 + sin(conj(q4))*conj(L4)*vq4^2 + sin(conj(q5))*conj(d5)*vq5^2 + aG3i - aq4*cos(conj(q4))*conj(L4) - aq3*cos(conj(q3))*conj(d3) - aq5*cos(conj(q5))*conj(d5);
    aG5j=cos(conj(q3))*conj(d3)*vq3^2 + cos(conj(q4))*conj(L4)*vq4^2 + cos(conj(q5))*conj(d5)*vq5^2 + aG3j + aq4*sin(conj(q4))*conj(L4) + aq3*sin(conj(q3))*conj(d3) + aq5*sin(conj(q5))*conj(d5);
    aG_mat(:,i)=[aG1i,aG1j,aG2i,aG2j,aG3i,aG3j,aG4i,aG4j,aG5i,aG5j]';
    
    % GRF
    aGi_list=aG_mat(1:2:end,i); aGj_list=aG_mat(2:2:end,i);
    GRF_mat(1,i)=1/2*sum(m_list.*aGi_list'); GRF_mat(3,i)=1/2*sum(m_list.*aGi_list');
    GRF_mat(2,i)=1/2*(sum(m_list.*aGj_list')+sum(m_list*g)); GRF_mat(4,i)=1/2*(sum(m_list.*aGj_list')+sum(m_list*g));
    
    % Dynamics:
    Pu_list=[P_mat(1:6,i)',P_mat(5:6,i)',P_mat(9:10,i)'];
    Pui_list=Pu_list(1:2:end); Puj_list=Pu_list(2:2:end); Gi_list=G_mat(1:2:end,i); Gj_list=G_mat(2:2:end,i); aGi_list=aG_mat(1:2:end,i); aGj_list=aG_mat(2:2:end,i);
    u_list=zeros(6,1);
    for k=1:5
        Pi=Pui_list(k); Pj=Puj_list(k);
        torqueG_list=-m_list*g.*(Gi_list-Pi);
        accel_list=m_list.*((Gi_list-Pi).*aGj_list-(Gj_list-Pj).*aGi_list)-I_list.*aq_list;
        GRF_P5=(P_mat(11,i)-Pi)*GRF_mat(4,i)-(P_mat(12,i)-Pj)*GRF_mat(3,i);
        u_list(k)=-GRF_P5-sum(torqueG_list(k:5))+sum(accel_list(k:5));
    end
    u_mat(:,i)=u_list;
end

var_list_guess((15*n+1):(27*n))=reshape(P_mat,1,12*n);
var_list_guess((27*n+1):(29*n))=reshape(vP0_mat,1,2*n);
var_list_guess((29*n+1):(31*n))=reshape(aP0_mat,1,2*n);
var_list_guess((31*n+1):(33*n))=reshape(vP5_mat,1,2*n);
var_list_guess((33*n+1):(35*n))=reshape(aP5_mat,1,2*n);
var_list_guess((35*n+1):(45*n))=reshape(G_mat,1,10*n);
var_list_guess((45*n+1):(55*n))=reshape(vG_mat,1,10*n);
var_list_guess((55*n+1):(65*n))=reshape(aG_mat,1,10*n);
var_list_guess((65*n+1):(69*n))=reshape(GRF_mat,1,4*n);
var_list_guess((69*n+1):(75*n))=reshape(u_mat,1,6*n);
%var_list_guess((69*n+1):(75*n))=sqrt(2e4/tf/6)*ones(1,6*n);
if iswalking
    var_list_guess([75*n+1,75*n+2])=[tf/2,tf/2]; % for walking
else
    var_list_guess([75*n+1,75*n+2])=[tf/2-sqrt((h_run/6+d_run/4)/2/g),tf/2+sqrt((h_run/6+d_run/4)/2/g)];
end


end
function [c,ceq] = running_nlcon(var_list,n,m_list,I_list,g,L_list,d_list,d_run,h_run,tf,h,ft_clr)
% Nonlinear conditions.
% P,G,vG,aG are symbolic expressions dependent on q,vq,aq.

ceq=[]; %ceq=0
c=[]; %c<=0

q_mat=reshape(var_list(1:(5*n)),5,n);
vq_mat=reshape(var_list((5*n+1):(10*n)),5,n);
aq_mat=reshape(var_list((10*n+1):(15*n)),5,n);
P_mat=reshape(var_list((15*n+1):(27*n)),12,n);
vP0_mat=reshape(var_list((27*n+1):(29*n)),2,n);
aP0_mat=reshape(var_list((29*n+1):(31*n)),2,n);
vP5_mat=reshape(var_list((31*n+1):(33*n)),2,n);
aP5_mat=reshape(var_list((33*n+1):(35*n)),2,n);
G_mat=reshape(var_list((35*n+1):(45*n)),10,n);
vG_mat=reshape(var_list((45*n+1):(55*n)),10,n);
aG_mat=reshape(var_list((55*n+1):(65*n)),10,n);
GRF_mat=reshape(var_list((65*n+1):(69*n)),4,n);
u_mat=reshape(var_list((69*n+1):(75*n)),6,n);

t1=var_list(75*n+1); t2=var_list(75*n+2);
n1=(n-1)/tf*t1+1; n2=(n-1)/tf*t2+1;

L1=L_list(1);L2=L_list(2);L3=L_list(3);L4=L_list(4);L5=L_list(5);
d1=d_list(1);d2=d_list(2);d3=d_list(3);d4=d_list(4);d5=d_list(5);

%kinematics and dynamics
for i=1:n % at every knot point
    
    q_list=q_mat(:,i); vq_list=vq_mat(:,i); aq_list=aq_mat(:,i); P_list=P_mat(:,i); vP0_list=vP0_mat(:,i); aP0_list=aP0_mat(:,i); vP5_list=vP5_mat(:,i); aP5_list=aP5_mat(:,i); G_list=G_mat(:,i); vG_list=vG_mat(:,i); aG_list=aG_mat(:,i); GRF_list=GRF_mat(:,i); u_list=u_mat(:,i);
    
    % Dynamics:
    Pu_list=[P_list(1:6)',P_list(5:6)',P_list(9:10)'];
    Pui_list=Pu_list(1:2:end); Puj_list=Pu_list(2:2:end); Gi_list=G_list(1:2:end); Gj_list=G_list(2:2:end); aGi_list=aG_list(1:2:end); aGj_list=aG_list(2:2:end); 
    for k=1:5
        Pi=Pui_list(k); Pj=Puj_list(k);
        torqueG_list=-m_list*g.*(Gi_list-Pi);
        accel_list=m_list.*((Gi_list-Pi).*aGj_list-(Gj_list-Pj).*aGi_list)-I_list.*aq_list;
        GRF_P5=(P_list(11)-Pi)*GRF_list(4)-(P_list(12)-Pj)*GRF_list(3);
        eq_dynamics=u_list(k)-u_list(6)+GRF_P5+sum(torqueG_list(k:5))-sum(accel_list(k:5));
        ceq=[ceq,eq_dynamics];
    end

    % Kinematics (body-based):
    % !!! PUT KINEMATICS EXPLICITLY HERE!!!
    G3i=G_list(5); G3j=G_list(6); vG3i=vG_list(5); vG3j=vG_list(6); aG3i=aG_list(5); aG3j=aG_list(6);
    
    q1=q_list(1); q2=q_list(2); q3=q_list(3); q4=q_list(4); q5=q_list(5);
    vq1=vq_list(1); vq2=vq_list(2); vq3=vq_list(3); vq4=vq_list(4); vq5=vq_list(5);
    aq1=aq_list(1); aq2=aq_list(2); aq3=aq_list(3); aq4=aq_list(4); aq5=aq_list(5);
    
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
    
    vP0i=vG3i - vq1*cos(conj(q1))*conj(L1) - vq2*cos(conj(q2))*conj(L2) - vq3*cos(conj(q3))*conj(d3);
    vP0j=vG3j + vq1*sin(conj(q1))*conj(L1) + vq2*sin(conj(q2))*conj(L2) + vq3*sin(conj(q3))*conj(d3);
    vP0=[vP0i,vP0j]';
    vP5i=vG3i - vq4*cos(conj(q4))*conj(L4) - vq5*cos(conj(q5))*conj(L5) - vq3*cos(conj(q3))*conj(d3);
    vP5j=vG3j + vq4*sin(conj(q4))*conj(L4) + vq5*sin(conj(q5))*conj(L5) + vq3*sin(conj(q3))*conj(d3);
    vP5=[vP5i,vP5j]';
    
    aP0i=sin(conj(q1))*conj(L1)*vq1^2 + sin(conj(q2))*conj(L2)*vq2^2 + sin(conj(q3))*conj(d3)*vq3^2 + aG3i - aq1*cos(conj(q1))*conj(L1) - aq2*cos(conj(q2))*conj(L2) - aq3*cos(conj(q3))*conj(d3);
    aP0j=cos(conj(q1))*conj(L1)*vq1^2 + cos(conj(q2))*conj(L2)*vq2^2 + cos(conj(q3))*conj(d3)*vq3^2 + aG3j + aq1*sin(conj(q1))*conj(L1) + aq2*sin(conj(q2))*conj(L2) + aq3*sin(conj(q3))*conj(d3);
    aP0=[aP0i,aP0j]';
    aP5i=sin(conj(q3))*conj(d3)*vq3^2 + sin(conj(q4))*conj(L4)*vq4^2 + sin(conj(q5))*conj(L5)*vq5^2 + aG3i - aq4*cos(conj(q4))*conj(L4) - aq5*cos(conj(q5))*conj(L5) - aq3*cos(conj(q3))*conj(d3);
    aP5j=cos(conj(q3))*conj(d3)*vq3^2 + cos(conj(q4))*conj(L4)*vq4^2 + cos(conj(q5))*conj(L5)*vq5^2 + aG3j + aq4*sin(conj(q4))*conj(L4) + aq5*sin(conj(q5))*conj(L5) + aq3*sin(conj(q3))*conj(d3);
    aP5=[aP5i,aP5j]';

    vG1i=vG3i - vq2*cos(conj(q2))*conj(L2) - vq1*cos(conj(q1))*conj(d1) - vq3*cos(conj(q3))*conj(d3);
    vG1j=vG3j + vq2*sin(conj(q2))*conj(L2) + vq1*sin(conj(q1))*conj(d1) + vq3*sin(conj(q3))*conj(d3);
    vG2i=vG3i - vq2*cos(conj(q2))*conj(d2) - vq3*cos(conj(q3))*conj(d3);
    vG2j=vG3j + vq2*sin(conj(q2))*conj(d2) + vq3*sin(conj(q3))*conj(d3);
    vG4i=vG3i - vq3*cos(conj(q3))*conj(d3) - vq4*cos(conj(q4))*conj(d4);
    vG4j=vG3j + vq3*sin(conj(q3))*conj(d3) + vq4*sin(conj(q4))*conj(d4);
    vG5i=vG3i - vq4*cos(conj(q4))*conj(L4) - vq3*cos(conj(q3))*conj(d3) - vq5*cos(conj(q5))*conj(d5);
    vG5j=vG3j + vq4*sin(conj(q4))*conj(L4) + vq3*sin(conj(q3))*conj(d3) + vq5*sin(conj(q5))*conj(d5);
    vG=[vG1i,vG1j,vG2i,vG2j,vG3i,vG3j,vG4i,vG4j,vG5i,vG5j]';
    
    aG1i=sin(conj(q1))*conj(d1)*vq1^2 + sin(conj(q2))*conj(L2)*vq2^2 + sin(conj(q3))*conj(d3)*vq3^2 + aG3i - aq2*cos(conj(q2))*conj(L2) - aq1*cos(conj(q1))*conj(d1) - aq3*cos(conj(q3))*conj(d3);
    aG1j=cos(conj(q1))*conj(d1)*vq1^2 + cos(conj(q2))*conj(L2)*vq2^2 + cos(conj(q3))*conj(d3)*vq3^2 + aG3j + aq2*sin(conj(q2))*conj(L2) + aq1*sin(conj(q1))*conj(d1) + aq3*sin(conj(q3))*conj(d3);
    aG2i=sin(conj(q2))*conj(d2)*vq2^2 + sin(conj(q3))*conj(d3)*vq3^2 + aG3i - aq2*cos(conj(q2))*conj(d2) - aq3*cos(conj(q3))*conj(d3);
    aG2j=cos(conj(q2))*conj(d2)*vq2^2 + cos(conj(q3))*conj(d3)*vq3^2 + aG3j + aq2*sin(conj(q2))*conj(d2) + aq3*sin(conj(q3))*conj(d3);
    aG4i=sin(conj(q3))*conj(d3)*vq3^2 + sin(conj(q4))*conj(d4)*vq4^2 + aG3i - aq3*cos(conj(q3))*conj(d3) - aq4*cos(conj(q4))*conj(d4);
    aG4j=cos(conj(q3))*conj(d3)*vq3^2 + cos(conj(q4))*conj(d4)*vq4^2 + aG3j + aq3*sin(conj(q3))*conj(d3) + aq4*sin(conj(q4))*conj(d4);
    aG5i=sin(conj(q3))*conj(d3)*vq3^2 + sin(conj(q4))*conj(L4)*vq4^2 + sin(conj(q5))*conj(d5)*vq5^2 + aG3i - aq4*cos(conj(q4))*conj(L4) - aq3*cos(conj(q3))*conj(d3) - aq5*cos(conj(q5))*conj(d5);
    aG5j=cos(conj(q3))*conj(d3)*vq3^2 + cos(conj(q4))*conj(L4)*vq4^2 + cos(conj(q5))*conj(d5)*vq5^2 + aG3j + aq4*sin(conj(q4))*conj(L4) + aq3*sin(conj(q3))*conj(d3) + aq5*sin(conj(q5))*conj(d5);
    aG=[aG1i,aG1j,aG2i,aG2j,aG3i,aG3j,aG4i,aG4j,aG5i,aG5j]';

    ceq=[ceq,(P_list-P)',(G_list-G)',(vG_list-vG)',(aG_list-aG)',(vP0_list-vP0)',(vP5_list-vP5)',(aP0_list-aP0)',(aP5_list-aP5)'];    
    
    % Terrain constraints: P0j>=h_jump*1(P0i-d_jump)
    P0i=P_list(1); P0j=P_list(2); P5i=P_list(11); P5j=P_list(12);
    h_terrain_0=h_run*heaviside(P0i-(d_run-0.1)+1e-6);
    h_terrain_1=h_run*heaviside(P1i-(d_run-0.1)+1e-6);
    h_terrain_4=h_run*heaviside(P4i-(d_run-0.1)+1e-6);
    h_terrain_5=h_run*heaviside(P5i-(d_run-0.1)+1e-6);
    c=[c,(h_terrain_0-P0j),(h_terrain_5-P5j),(h_terrain_1-P1j),(h_terrain_4-P4j)];
        
    % Phase constraints:
    P0i=P_list(1); P0j=P_list(2); vP0i=vP0_list(1); vP0j=vP0_list(2); aP0i=aP0_list(1); aP0j=aP0_list(2);
    P5i=P_list(11); P5j=P_list(12); vP5i=vP5_list(1); vP5j=vP5_list(2); aP5i=aP5_list(1); aP5j=aP5_list(2);
    if i<=n1
        % Contact phase1 constraints:
        ceq=[ceq,GRF_list(3:4)',u_list(6),P0i,P0j,vP0i,vP0j,aP0i,aP0j];
        %c=[c,0,h_terrain_5+ft_clr-P5j];
    elseif (i>n1)&&(i<=n2)
        % Aerial Phase constraints
        ceq=[ceq,GRF_list(1:4)',u_list(1),u_list(6),0,0,0];
        %c=[c,h_terrain_0+ft_clr-P0j,h_terrain_5+ft_clr-P5j];
    elseif i>n2
        % Contact phase2 constraints:
        ceq=[ceq,GRF_list(1:2)',u_list(1),P5i-d_run,P5j-h_run,vP5i,vP5j,aP5i,aP5j];
        %c=[c,h_terrain_0+ft_clr-P0j,0];
    end
    
    % Collocation on body G3 for aerial phase and phase transition
    if i<=(n-1)
        ceq_co1=G_mat(5,i+1)-G_mat(5,i)-h/2*(vG_mat(5,i)+vG_mat(5,i+1)); % G3i and vG3i
        ceq_co2=G_mat(6,i+1)-G_mat(6,i)-h/2*(vG_mat(6,i)+vG_mat(6,i+1)); % G3j and vG3j
        ceq_co3=vG_mat(5,i+1)-vG_mat(5,i)-h/2*(aG_mat(5,i)+aG_mat(5,i+1)); % vG3i and aG3i
        ceq_co4=vG_mat(6,i+1)-vG_mat(6,i)-h/2*(aG_mat(6,i)+aG_mat(6,i+1)); % vG3j and aG3j
        if (i>(n1-1))&&(i<=n2)
            ceq=[ceq,ceq_co1,ceq_co2,ceq_co3,ceq_co4];
        else
            ceq=[ceq,0,0,0,0];
        end
    end

end

end
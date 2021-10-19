function [c,ceq] = standing_jump_nlcon_contact_implicit(var_list,n,m_list,I_list,g,L_list,d_list,d_jump,h_jump,miu,tf,h)
% Nonlinear conditions.
% P,G,vG,aG are symbolic expressions dependent on q,vq,aq.

ceq=[]; %ceq=0
c=[]; %c<=0

q_mat=reshape(var_list(1:(3*n)),3,n);
vq_mat=reshape(var_list((3*n+1):(6*n)),3,n);
aq_mat=reshape(var_list((6*n+1):(9*n)),3,n);
P_mat=reshape(var_list((9*n+1):(17*n)),8,n);
vP0_mat=reshape(var_list((17*n+1):(19*n)),2,n);
aP0_mat=reshape(var_list((19*n+1):(21*n)),2,n);
G_mat=reshape(var_list((21*n+1):(27*n)),6,n);
vG_mat=reshape(var_list((27*n+1):(33*n)),6,n);
aG_mat=reshape(var_list((33*n+1):(39*n)),6,n);
u_mat=reshape(var_list((39*n+1):(42*n)),3,n);

L1=L_list(1);L2=L_list(2);L3=L_list(3);
d1=d_list(1);d2=d_list(2);d3=d_list(3);

%standing jump kinematics and dynamics
for i=1:n % at every knot point
    
    q_list=q_mat(:,i); vq_list=vq_mat(:,i); aq_list=aq_mat(:,i); P_list=P_mat(:,i); vP0_list=vP0_mat(:,i); aP0_list=aP0_mat(:,i); G_list=G_mat(:,i); vG_list=vG_mat(:,i); aG_list=aG_mat(:,i); u_list=u_mat(:,i);
    
    % Dynamics:
    Pi_list=P_list(1:2:end); Pj_list=P_list(2:2:end); Gi_list=G_list(1:2:end); Gj_list=G_list(2:2:end); aGi_list=aG_list(1:2:end); aGj_list=aG_list(2:2:end); 
    for k=1:length(u_list)
        Pi=Pi_list(k); Pj=Pj_list(k);
        torqueG_list=-m_list*g.*(Gi_list-Pi);
        accel_list=m_list.*((Gi_list-Pi).*aGj_list-(Gj_list-Pj).*aGi_list)-I_list.*aq_list;
        eq_dynamics=u_list(k)+sum(torqueG_list(k:3))-sum(accel_list(k:3));
        ceq=[ceq,eq_dynamics];
    end

    % Kinematics:
    % !!! PUT KINEMATICS EXPLICITLY HERE!!!
    P0i=P_list(1); P0j=P_list(2); vP0i=vP0_mat(1,i); vP0j=vP0_mat(2,i); aP0i=aP0_mat(1,i); aP0j=aP0_mat(2,i);
    
    q1=q_list(1); q2=q_list(2); q3=q_list(3);
    vq1=vq_list(1); vq2=vq_list(2); vq3=vq_list(3);
    aq1=aq_list(1); aq2=aq_list(2); aq3=aq_list(3);
    
    P1i=P0i+L1*sin(q1); P1j=P0j+L1*cos(q1);
    P2i=P1i+L2*sin(q2); P2j=P1j+L2*cos(q2);
    P3i=P2i+L3*sin(q3); P3j=P2j+L3*cos(q3);
    P=[P0i,P0j,P1i,P1j,P2i,P2j,P3i,P3j]';

    G1i=P0i+d1*sin(q1); G1j=P0j+d1*cos(q1);
    G2i=P1i+d2*sin(q2); G2j=P1j+d2*cos(q2);
    G3i=P2i+d3*sin(q3); G3j=P2j+d3*cos(q3);
    G=[G1i,G1j,G2i,G2j,G3i,G3j]';

    vG1i=vP0i + vq1*cos(conj(q1))*conj(d1);
    vG1j=vP0j - vq1*sin(conj(q1))*conj(d1);
    vG2i=vP0i + vq1*cos(conj(q1))*conj(L1) + vq2*cos(conj(q2))*conj(d2);
    vG2j=vP0j - vq1*sin(conj(q1))*conj(L1) - vq2*sin(conj(q2))*conj(d2);
    vG3i=vP0i + vq1*cos(conj(q1))*conj(L1) + vq2*cos(conj(q2))*conj(L2) + vq3*cos(conj(q3))*conj(d3);
    vG3j=vP0j - vq1*sin(conj(q1))*conj(L1) - vq2*sin(conj(q2))*conj(L2) - vq3*sin(conj(q3))*conj(d3);
    vG=[vG1i,vG1j,vG2i,vG2j,vG3i,vG3j]';
    
    aG1i=- sin(conj(q1))*conj(d1)*vq1^2 + aP0i + aq1*cos(conj(q1))*conj(d1);
    aG1j=- cos(conj(q1))*conj(d1)*vq1^2 + aP0j - aq1*sin(conj(q1))*conj(d1);
    aG2i=- sin(conj(q1))*conj(L1)*vq1^2 - sin(conj(q2))*conj(d2)*vq2^2 + aP0i + aq1*cos(conj(q1))*conj(L1) + aq2*cos(conj(q2))*conj(d2);
    aG2j=- cos(conj(q1))*conj(L1)*vq1^2 - cos(conj(q2))*conj(d2)*vq2^2 + aP0j - aq1*sin(conj(q1))*conj(L1) - aq2*sin(conj(q2))*conj(d2);
    aG3i=- sin(conj(q1))*conj(L1)*vq1^2 - sin(conj(q2))*conj(L2)*vq2^2 - sin(conj(q3))*conj(d3)*vq3^2 + aP0i + aq1*cos(conj(q1))*conj(L1) + aq2*cos(conj(q2))*conj(L2) + aq3*cos(conj(q3))*conj(d3);
    aG3j=- cos(conj(q1))*conj(L1)*vq1^2 - cos(conj(q2))*conj(L2)*vq2^2 - cos(conj(q3))*conj(d3)*vq3^2 + aP0j - aq1*sin(conj(q1))*conj(L1) - aq2*sin(conj(q2))*conj(L2) - aq3*sin(conj(q3))*conj(d3);
    aG=[aG1i,aG1j,aG2i,aG2j,aG3i,aG3j]';

    ceq=[ceq,(P_list-P)',(G_list-G)',(vG_list-vG)',(aG_list-aG)'];
    
    % Terrain constraints: P0j>h_jump*1(P0i-d_jump)
    if d_jump>=0
        h_terrain=h_jump*heaviside(P0i-(d_jump-0.1)+1e-3);
    elseif d_jump<0
        h_terrain=h_jump*(1-heaviside(P0i-(d_jump-0.1)+1e-3));
    end
    c=[c,(h_terrain-P0j)];
    
    % Phase constraints:
    friction=sum(m_list.*aGi_list');
    normal=sum(m_list.*aGj_list')+sum(m_list*g);
    P0i=P_list(1); P0j=P_list(2); vP0i=vP0_list(1); vP0j=vP0_list(2); aP0i=aP0_list(1); aP0j=aP0_list(2);
    if P0j<=(h_terrain+1e-3)
        % Contact phase constraints:
        ceq=[ceq,vP0i,vP0j,aP0i,aP0j];
        c=[c,(1e-6-normal),(abs(friction)-miu*normal)];
    elseif P0j>(h_terrain+1e-3)
        % Aerial Phase constraints
        ceq=[ceq,friction,normal,u_list(1),0];
        c=[c,0,0];
    end
    
    % Collocation on body G3 for aerial phase and phase transition
    if i<=(n-1)
        ceq_co1=G_mat(5,i+1)-G_mat(5,i)-h/2*(vG_mat(5,i)+vG_mat(5,i+1)); % G3i and vG3i
        ceq_co2=G_mat(6,i+1)-G_mat(6,i)-h/2*(vG_mat(6,i)+vG_mat(6,i+1)); % G3j and vG3j
        ceq_co3=vG_mat(5,i+1)-vG_mat(5,i)-h/2*(aG_mat(5,i)+aG_mat(5,i+1)); % vG3i and aG3i
        ceq_co4=vG_mat(6,i+1)-vG_mat(6,i)-h/2*(aG_mat(6,i)+aG_mat(6,i+1)); % vG3j and aG3j
        if (P_mat(2,i)>(h_terrain+1e-3))||(P_mat(2,i+1)>(h_terrain+1e-3))
            ceq=[ceq,ceq_co1,ceq_co2,ceq_co3,ceq_co4];
        else
            ceq=[ceq,0,0,0,0];
        end
    end

end

end
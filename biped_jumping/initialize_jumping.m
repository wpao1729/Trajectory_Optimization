function var_list_guess=initialize_jumping(n,qi_list,qf_list,L_list,d_list,m_list,I_list,g,d_jump,h_jump,t_list,isflip)

var_list_guess=zeros(1,42*n+2);

L1=L_list(1);L2=L_list(2);L3=L_list(3);
d1=d_list(1);d2=d_list(2);d3=d_list(3);

q_mat=zeros(3,n); vq_mat=zeros(3,n); aq_mat=zeros(3,n);
tf=t_list(end);
for k=1:3
    a=-6*(qf_list(k)-qi_list(k))/(tf^3); 
    q_mat(k,:)=a/3*(t_list.^3)-a*tf/2*(t_list.^2)+qi_list(k);
    vq_mat(k,:)=a*(t_list.^2)-a*tf*t_list;
    aq_mat(k,:)=2*a*t_list-a*tf;
end

P0_mat=zeros(2,n); vP0_mat=zeros(2,n); aP0_mat=zeros(2,n);
P0_mat(1,:)=linspace(0,d_jump,n);
% Calculate a parabola for P0j(t) initial guess
if isflip
    dh=0.75;
else
    dh=0.3;
end
b=2*(h_jump+dh)*(1+sqrt(dh/(h_jump+dh)))/tf;
a=-b^2/(4*(h_jump+dh));
P0_mat(2,:)=a*t_list.^2+b*t_list;
vP0_mat(1,:)=d_jump/tf;
vP0_mat(2,:)=2*a*t_list+b;
% aP0i=0;
aP0_mat(2,:)=2*a;

var_list_guess(1:(3*n))=reshape(q_mat,1,3*n);
var_list_guess((3*n+1):(6*n))=reshape(vq_mat,1,3*n);
var_list_guess((6*n+1):(9*n))=reshape(aq_mat,1,3*n);
var_list_guess((17*n+1):(19*n))=reshape(vP0_mat,1,2*n);
var_list_guess((19*n+1):(21*n))=reshape(aP0_mat,1,2*n);

P_mat=zeros(8,n); G_mat=zeros(6,n); vG_mat=zeros(6,n); aG_mat=zeros(6,n); u_mat=zeros(3,n);
for i=1:n % at every knot point

    q_list=q_mat(:,i); vq_list=vq_mat(:,i); aq_list=aq_mat(:,i);
    P0i=P0_mat(1,i); P0j=P0_mat(2,i);

    % Kinematics:
    vP0i=vP0_mat(1,i); vP0j=vP0_mat(2,i); aP0i=aP0_mat(1,i); aP0j=aP0_mat(2,i);
    
    q1=q_list(1); q2=q_list(2); q3=q_list(3);
    vq1=vq_list(1); vq2=vq_list(2); vq3=vq_list(3);
    aq1=aq_list(1); aq2=aq_list(2); aq3=aq_list(3);
    
    P1i=P0i+L1*sin(q1); P1j=P0j+L1*cos(q1);
    P2i=P1i+L2*sin(q2); P2j=P1j+L2*cos(q2);
    P3i=P2i+L3*sin(q3); P3j=P2j+L3*cos(q3);
    P_mat(:,i)=[P0i,P0j,P1i,P1j,P2i,P2j,P3i,P3j]';

    G1i=P0i+d1*sin(q1); G1j=P0j+d1*cos(q1);
    G2i=P1i+d2*sin(q2); G2j=P1j+d2*cos(q2);
    G3i=P2i+d3*sin(q3); G3j=P2j+d3*cos(q3);
    G_mat(:,i)=[G1i,G1j,G2i,G2j,G3i,G3j]';

    vG1i=vP0i + vq1*cos(conj(q1))*conj(d1);
    vG1j=vP0j - vq1*sin(conj(q1))*conj(d1);
    vG2i=vP0i + vq1*cos(conj(q1))*conj(L1) + vq2*cos(conj(q2))*conj(d2);
    vG2j=vP0j - vq1*sin(conj(q1))*conj(L1) - vq2*sin(conj(q2))*conj(d2);
    vG3i=vP0i + vq1*cos(conj(q1))*conj(L1) + vq2*cos(conj(q2))*conj(L2) + vq3*cos(conj(q3))*conj(d3);
    vG3j=vP0j - vq1*sin(conj(q1))*conj(L1) - vq2*sin(conj(q2))*conj(L2) - vq3*sin(conj(q3))*conj(d3);
    vG_mat(:,i)=[vG1i,vG1j,vG2i,vG2j,vG3i,vG3j]';
    
    aG1i=- sin(conj(q1))*conj(d1)*vq1^2 + aP0i + aq1*cos(conj(q1))*conj(d1);
    aG1j=- cos(conj(q1))*conj(d1)*vq1^2 + aP0j - aq1*sin(conj(q1))*conj(d1);
    aG2i=- sin(conj(q1))*conj(L1)*vq1^2 - sin(conj(q2))*conj(d2)*vq2^2 + aP0i + aq1*cos(conj(q1))*conj(L1) + aq2*cos(conj(q2))*conj(d2);
    aG2j=- cos(conj(q1))*conj(L1)*vq1^2 - cos(conj(q2))*conj(d2)*vq2^2 + aP0j - aq1*sin(conj(q1))*conj(L1) - aq2*sin(conj(q2))*conj(d2);
    aG3i=- sin(conj(q1))*conj(L1)*vq1^2 - sin(conj(q2))*conj(L2)*vq2^2 - sin(conj(q3))*conj(d3)*vq3^2 + aP0i + aq1*cos(conj(q1))*conj(L1) + aq2*cos(conj(q2))*conj(L2) + aq3*cos(conj(q3))*conj(d3);
    aG3j=- cos(conj(q1))*conj(L1)*vq1^2 - cos(conj(q2))*conj(L2)*vq2^2 - cos(conj(q3))*conj(d3)*vq3^2 + aP0j - aq1*sin(conj(q1))*conj(L1) - aq2*sin(conj(q2))*conj(L2) - aq3*sin(conj(q3))*conj(d3);
    aG_mat(:,i)=[aG1i,aG1j,aG2i,aG2j,aG3i,aG3j]';
    
    % Dynamics:
    Pi_list=P_mat(1:2:end,i); Pj_list=P_mat(2:2:end,i); Gi_list=G_mat(1:2:end,i); Gj_list=G_mat(2:2:end,i); aGi_list=aG_mat(1:2:end,i); aGj_list=aG_mat(2:2:end,i); 
    u_list=zeros(3,1);
    for k=1:length(u_list)
        Pi=Pi_list(k); Pj=Pj_list(k);
        torqueG_list=-m_list*g.*(Gi_list-Pi);
        accel_list=m_list.*((Gi_list-Pi).*aGj_list-(Gj_list-Pj).*aGi_list)-I_list.*aq_list;
        u_list(k)=-sum(torqueG_list(k:3))+sum(accel_list(k:3));
    end
    u_mat(:,i)=u_list;
end

var_list_guess((9*n+1):(17*n))=reshape(P_mat,1,8*n);
var_list_guess((21*n+1):(27*n))=reshape(G_mat,1,6*n);
var_list_guess((27*n+1):(33*n))=reshape(vG_mat,1,6*n);
var_list_guess((33*n+1):(39*n))=reshape(aG_mat,1,6*n);

% Ensure the cost of initial guess is high enough for better performance of subproblem algorithm cg
h=tf/(n-1);
effort_list=sum(u_mat.^2);
cost=h/2*(effort_list(1)+effort_list(end)+2*sum(effort_list(2:(end-1))));

if isflip
    if cost>=2e5
        var_list_guess((39*n+1):(42*n))=reshape(u_mat,1,3*n);
    else
        var_list_guess((39*n+1):(42*n))=reshape(u_mat*sqrt(3e5/cost),1,3*n);
    end
    var_list_guess([42*n+1,42*n+2])=[tf/2-sqrt((h_jump+2.5)/2/g),tf/2+sqrt((h_jump+2.5)/2/g)];
else
    if cost>=5e4
        var_list_guess((39*n+1):(42*n))=reshape(u_mat,1,3*n);
    else
        var_list_guess((39*n+1):(42*n))=sqrt(8e4/tf/3)*ones(1,3*n);
    end
    var_list_guess([42*n+1,42*n+2])=[tf/2-sqrt((h_jump+0.5)/2/g),tf/2+sqrt((h_jump+0.5)/2/g)]; % typically +0.2 ~ +0.5
end

end
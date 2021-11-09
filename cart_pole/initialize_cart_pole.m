function var_list_guess = initialize_cart_pole (n,x1,xf,theta1,thetaf,m_cart,m_pole,L)

T_guess=2; 
t_list_guess=linspace(0,T_guess,n);
a_x=-6*(xf-x1)/(T_guess^3);
x_list_guess=a_x/3*(t_list_guess.^3)-a_x*T_guess/2*(t_list_guess.^2)+x1;
v_list_guess=a_x*(t_list_guess.^2)-a_x*T_guess*t_list_guess;
a_list_guess=2*a_x*t_list_guess-a_x*T_guess;

a_theta=-6*(thetaf-theta1)/(T_guess^3);
theta_list_guess=a_theta/3*(t_list_guess.^3)-a_theta*T_guess/2*(t_list_guess.^2)+theta1;
omega_list_guess=a_theta*(t_list_guess.^2)-a_theta*T_guess*t_list_guess;
alpha_list_guess=2*a_theta*t_list_guess-a_theta*T_guess;

u_list_guess=(m_cart+m_pole)*a_list_guess-m_pole*L*sin(theta_list_guess).*alpha_list_guess-m_pole*L*cos(theta_list_guess).*(omega_list_guess.^2);
var_list_guess=[x_list_guess,v_list_guess,a_list_guess,theta_list_guess,omega_list_guess,alpha_list_guess,u_list_guess,T_guess]; % Initial guess of the solution

end
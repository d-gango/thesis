clear all
% parameters
global m g F_ext  R t_ss a bb c V_star mu_star
m = 10;
g = 9.81;
F_ext = 5;
a = 0.0349;
bb = 0.0489;
c = 10e-3;
R = 100;
V_star = 10e-6;
mu_star = 0.369;
t_ss = 100;
% initial conditions
x0 = [10e-6;0;0;];
T = 20;
% solve ODE
[t,x] = ode15s(@eom, [0,T], x0);
% calculate friction coeff
mu = a.*asinh((x(:,2)./V_star)./2.*exp(mu_star+bb.*log(c+x(:,3))./a));
% friction force
Ff = mu*m*g;
Fsum = F_ext - Ff;


figure
plot(t,x(:,1), 'o-')
title('position')
figure
plot(t,x(:,2), 'o-')
title('velocity')
figure
plot(t,x(:,3), 'o-')
title('$\bar{\phi}$ ("surface roughness")', 'Interpreter', 'latex')
figure
plot(t,mu, 'o-')
title('mu')
figure
plot(t,Ff, 'o-')
title('friction force')
figure
plot(t,Fsum, 'o-')
title('net force')

function xdot = eom(t,x)
% x = [pos; vel; phibar]
global m F_ext a bb c t_ss R V_star mu_star g
phi_bar = x(3);
Fn = m*g;
% friction force
v = x(2) / V_star;
mu = a*asinh(v/2*exp(mu_star+bb*log(c+phi_bar)/a));
Ff = mu*Fn;
% evolution of contact surface roughness
phi_bardot = -(1+R*v)/t_ss * sinh((R*v*phi_bar-(1+R-phi_bar))/(1+R*v));
Fsum = F_ext - Ff;
vel = x(2);
acc = Fsum/m;
xdot = [vel;acc;phi_bardot];
end

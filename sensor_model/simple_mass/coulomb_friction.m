clear all
% parameters
global m k b mu v_f g
m = 1;
k = 1;
b = 1;
g = 10;
mu = 0.5;
v_f = @(t) heaviside(t);
% initial conditions
x0 = [0;0];
T = 15;
% solve ODE
[t,x] = ode15s(@eom, [0,T], x0);
% calculate friction force
v_rel = x(:,2) - v_f(t); %velocity of the body relative to the ground
Ff = -mu*tanh(10*v_rel);

figure
plot(t,x(:,1), 'o-')
title('position')
figure
plot(t,x(:,2), 'o-')
title('velocity')
figure
plot(t,v_rel, 'o-')
title('relative velocity')
figure
plot(t,Ff, 'o-')
title('friction force')

function xdot = eom(t,x)
% x = [pos; vel]
global m k b mu v_f g
Fk = -k*x(1);
Fb = -b*x(2);
vrel = x(2) - v_f(t);
Ff = -mu*tanh(10*vrel)*m*g;
Fsum = Fk+Fb+Ff;

vel = x(2);
acc = Fsum/m;
xdot = [vel;acc];
end

clear all
% parameters
global m1 m2 k b mu F_ext
m1 = 10;
m2 = 1;
k = 10;
b = 10;
mu = 0.5;
F_ext = 0.5;
% initial conditions
x0 = [0;0;0;0];
T = 15;
% solve ODE
[t,x] = ode15s(@eom, [0,T], x0);
% calculate friction force
v_rel = x(:,2) - x(:,4); %velocity of body 1 relative to body 2
Ff = -mu*tanh(10*v_rel);

figure
plot(t,x(:,1), 'o-')
title('position 1')
figure
plot(t,x(:,2), 'o-')
title('velocity 1')
figure
plot(t,x(:,3), 'o-')
title('position 2')
figure
plot(t,x(:,4), 'o-')
title('velocity 2')
figure
plot(t,v_rel, 'o-')
title('relative velocity')
figure
plot(t,Ff, 'o-')
title('friction force')

function xdot = eom(t,x)
% x = [pos; vel]
global m1 m2 k b mu F_ext
Fk = -k*x(1);
Fb = -b*x(2);
vrel = x(2) - x(4);
Ff = -mu*tanh(10*vrel);
Fsum1 = Fk+Fb+Ff;
Fsum2 = F_ext - Ff;

vel1 = x(2);
vel2 = x(4);
acc1 = Fsum1/m1;
acc2 = Fsum2/m2;
xdot = [vel1;acc1;vel2;acc2];
end

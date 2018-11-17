clear all
close all

compare = 0; % 1 - use ode15s and ode45, 0 - ode15s only

% parameters
global m k b mu v_f g smooth_friction
smooth_friction = @(vrel) tanh(100*vrel);
m = 10;
k = 100;
b = 10;
g = 9.81;
mu = 0.5134;
v_f = @(t) 1;
% initial conditions
x0 = [0;0];
T = 10;

switch compare
    case 0
        % solve ODE
        options = odeset('RelTol',1e-10,'AbsTol',1e-12,'InitialStep',1e-20);
        [t,x] = ode15s(@eom, [0,T], x0, options);
        % calculate friction force
        v_rel = v_f(t) - x(:,2); %velocity of the body relative to the ground
        Ff = m*g*mu*smooth_friction(v_rel);
        % net force
        Fk = -k*x(:,1);
        Fb = -b*x(:,2);
        Fsum = Fk+Fb+Ff;
        
        figure
        plot(t,x(:,1), '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('position')
        ylabel('x [m]'); xlabel('t [s]');
        
        figure
        plot(t,x(:,2), '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('velocity')
        ylabel('v [m/s]'); xlabel('t [s]');
        
        figure
        plot(t,v_rel, '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('relative velocity')
        ylabel('v_{rel} [m/s]'); xlabel('t [s]');
        
        figure
        plot(t,Ff, '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('friction force')
        ylabel('F_f [N]'); xlabel('t [s]');
        
        figure
        plot(t,Fsum, '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('net force')
        ylabel('F_s [N]'); xlabel('t [s]');

    case 1
        % solve ODE
        [t,x] = ode15s(@eom, [0,T], x0);
        [t45,x45] = ode45(@eom, [0,T], x0);
        % calculate friction force
        v_rel = v_f(t) - x(:,2); %velocity of the body relative to the ground
        v_rel45 = v_f(t45) - x45(:,2); %velocity of the body relative to the ground
        Ff = -m*g*mu*smooth_friction(v_rel);
        Ff45 = -m*g*mu*smooth_friction(v_rel45);
        
        figure
        plot(t,x(:,1), '.-', t45,x45(:,1), 'r.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('position')
        ylabel('x [m]'); xlabel('t [s]');
        legend('ode15s', 'ode45')
        
        figure
        plot(t,x(:,2), '.-', t45,x45(:,2), 'r.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('velocity')
        ylabel('v [m/s]'); xlabel('t [s]');
        legend('ode15s', 'ode45')
        
        figure
        plot(t,v_rel, '.-', t45,v_rel45, 'r.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('relative velocity')
        ylabel('v_{rel} [m/s]'); xlabel('t [s]');
        legend('ode15s', 'ode45')
        
        figure
        plot(t,Ff, '.-', t45,Ff45, 'r.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('friction force')
        ylabel('F_f [N]'); xlabel('t [s]');
        legend('ode15s', 'ode45')
end

coulomb.x = x;
coulomb.Ff = Ff;
coulomb.Fk = Fk;
coulomb.Fb = Fb;
coulomb.Fsum = Fsum;
coulomb.mu = mu;
coulomb.t = t;
coulomb.vf = v_f;
coulomb.v_rel = v_rel;
save coulomb_data_v.mat coulomb


function xdot = eom(t,x)
% x = [pos; vel]
global m k b mu v_f g smooth_friction
Fk = -k*x(1);
Fb = -b*x(2);
vrel = v_f(t) - x(2);
Ff = mu*smooth_friction(vrel)*m*g;
Fsum = Fk+Fb+Ff;

vel = x(2);
acc = Fsum/m;
xdot = [vel;acc];
end
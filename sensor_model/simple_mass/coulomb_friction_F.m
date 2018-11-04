clear all

compare = 0; % 1 - use ode15s and ode45, 0 - ode15s only

% parameters
global m mu F_ext g smooth_friction
smooth_friction = @(vrel) tanh(100*vrel);
m = 10;
g = 9.81;
mu = 0.4047;
F_ext = @(t) 10*t;
% initial conditions
x0 = [0;0];
T = 10;

switch compare
    case 0
        % solve ODE
        options = odeset('RelTol',1e-10,'AbsTol',1e-12,'InitialStep',1e-20);
        [t,x] = ode15s(@eom, [0,T], x0, options);
        % calculate friction force
        Ff = mu*smooth_friction(x(:,2))*m*g;
        Fsum = F_ext(t) - Ff;
        
        figure
        plot(t,x(:,1), '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('position')
        ylabel('x [m]'); xlabel('t [s]');
        
        figure
        plot(t,x(:,2), '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('velocity')
        ylabel('v [m/s]'); xlabel('t [s]');
        
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
        % calculate friction force
        Ff = mu*smooth_friction(x(:,2))*m*g;
        Fsum = F_ext(t) - Ff;
        
        % solve ODE45
        [t45,x45] = ode45(@eom, [0,T], x0);
        % calculate friction force
        Ff45 = mu*smooth_friction(x45(:,2))*m*g;
        Fsum45 = F_ext(t45) - Ff45;
        
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
        plot(t,Ff, '.-', t45,Ff45, 'r.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('friction force')
        ylabel('F_f [N]'); xlabel('t [s]');
        legend('ode15s', 'ode45')
        
        figure
        plot(t,Fsum, '.-', t45,Fsum45, 'r.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('net force')
        ylabel('F_s [N]'); xlabel('t [s]');
        legend('ode15s', 'ode45')
end

coulomb.x = x;
coulomb.Ff = Ff;
coulomb.Fsum = Fsum;
coulomb.mu = mu;
coulomb.t = t;
coulomb.Fext = F_ext;
save coulomb_data_F.mat coulomb


function xdot = eom(t,x)
% x = [pos; vel]
global m mu F_ext g smooth_friction
vel = x(2);
Fn = m*g;
Ff = -mu*smooth_friction(vel)*Fn;
Fsum = F_ext(t) + Ff;

acc = Fsum/m;
xdot = [vel;acc];
end

clear all

compare = 0; % 1 - use ode15s and ode45, 0 - ode15s only

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
x0 = [0;0;101;];
T = 10;

switch compare
    case 0
        % solve ODE
        [t,x] = ode15s(@eom, [0,T], x0);
        % calculate friction coeff
        mu = a.*asinh((x(:,2)./V_star)./2.*exp(mu_star+bb.*log(c+x(:,3))./a));
        % friction force
        Ff = mu*m*g;
        Fsum = F_ext - Ff;
        
        figure
        plot(t,x(:,1), '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('position')
        ylabel('x [m]'); xlabel('t [s]');
        
        figure
        plot(t,x(:,2), '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('velocity')
        ylabel('v [m/s]'); xlabel('t [s]');
        
        figure
        plot(t,x(:,3), '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('"surface roughness"')
        ylabel('$\bar{\phi}$ [1]', 'Interpreter', 'latex'); xlabel('t [s]', 'Interpreter', 'latex');
        
        figure
        plot(t,mu, '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('friction coeffitient')
        ylabel('\mu [1]'); xlabel('t [s]');
        
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
        % calculate friction coeff
        mu = a.*asinh((x(:,2)./V_star)./2.*exp(mu_star+bb.*log(c+x(:,3))./a));
        % friction force
        Ff = mu*m*g;
        Fsum = F_ext - Ff;
        
        % solve ODE45
        options = odeset('RelTol',1e-10,'AbsTol',1e-12, 'InitialStep',1e-20);
        [t45,x45] = ode45(@eom, [0,T], x0, options);
        % calculate friction coeff
        mu45 = a.*asinh((x45(:,2)./V_star)./2.*exp(mu_star+bb.*log(c+x45(:,3))./a));
        % friction force
        Ff45 = mu45*m*g;
        Fsum45 = F_ext - Ff45;
        
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
        plot(t,x(:,3), '.-', t45,x45(:,3), 'r.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('"surface roughness"')
        ylabel('$\bar{\phi}$ [1]', 'Interpreter', 'latex'); xlabel('t [s]', 'Interpreter', 'latex');
        legend('ode15s', 'ode45')
        
        figure
        plot(t,mu, '.-', t45,mu45, 'r.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('friction coeffitient')
        ylabel('\mu [1]'); xlabel('t [s]');
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

function xdot = eom(t,x)
t
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

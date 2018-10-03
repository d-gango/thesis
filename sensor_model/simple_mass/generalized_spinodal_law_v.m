clear all

compare = 0; % 1 - use ode15s and ode45, 0 - ode15s only

% parameters
global m k b v_f g a bb c R V_star mu_star t_ss
m = 10;
k = 10;
b = 10;
g = 9.81;
v_f = @(t) 1;
a = 0.0349;
bb = 0.0489;
c = 10e-3;
R = 100;
V_star = 10e-6;
mu_star = 0.369;
t_ss = 100;

% initial conditions
x0 = [0;0;101];
T = 20;

switch compare
    case 0
        % solve ODE
        [t,x] = ode15s(@eom, [0,T], x0);
        % calculate friction force
        v_rel = v_f(t) - x(:,2); %velocity of the ground relative to the body
        mu = a.*asinh((v_rel./V_star)./2.*exp(mu_star+bb.*log(c+x(:,3))./a));
        % friction force
        Ff = mu*m*g;
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
        plot(t,x(:,3), '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('"surface roughness"')
        ylabel('$\bar{\phi}$ [1]', 'Interpreter', 'latex'); xlabel('t [s]', 'Interpreter', 'latex');
        
        figure
        plot(t,mu, '.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('friction coeffitient')
        ylabel('\mu [1]'); xlabel('t [s]');
        
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
        % calculate friction force
        v_rel = v_f(t) - x(:,2); %velocity of the ground relative to the body
        mu = a.*asinh((v_rel./V_star)./2.*exp(mu_star+bb.*log(c+x(:,3))./a));
        % friction force
        Ff = mu*m*g;
        % net force
        Fk = -k*x(:,1);
        Fb = -b*x(:,2);
        Fsum = Fk+Fb+Ff;
        
        % solve ODE45
        [t45,x45] = ode45(@eom, [0,T], x0);
        % calculate friction force
        v_rel45 = v_f(t) - x45(:,2); %velocity of the ground relative to the body
        mu45 = a.*asinh((v_rel45./V_star)./2.*exp(mu_star+bb.*log(c+x45(:,3))./a));
        % friction force
        Ff45 = mu45*m*g;
        % net force
        Fk45 = -k*x45(:,1);
        Fb45 = -b*x45(:,2);
        Fsum45 = Fk45+Fb45+Ff45;
        
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
        plot(t,v_rel, '.-', t45,v_rel45, 'r.-', 'LineWidth', 2, 'MarkerSize', 15)
        title('relative velocity')
        ylabel('v_{rel} [m/s]'); xlabel('t [s]');
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
% x = [pos; vel; phi_bar]
global m k b v_f g a bb c R V_star mu_star t_ss
phi_bar = x(3);

Fk = -k*x(1);
Fb = -b*x(2);
Fn = m*g;
% friction force
vrel = v_f(t) - x(2); % floor velocity relative to the body
v = vrel / V_star;
mu = a*asinh(v/2*exp(mu_star+bb*log(c+phi_bar)/a));
Ff = mu*Fn;
% evolution of contact surface roughness
phi_bardot = -(1+R*v)/t_ss * sinh((R*v*phi_bar-(1+R-phi_bar))/(1+R*v));

Fsum = Fk+Fb+Ff;
vel = x(2);
acc = Fsum/m;
xdot = [vel;acc;phi_bardot];
end

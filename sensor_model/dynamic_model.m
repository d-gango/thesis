%% equation of motion
tic
clear all
par = param(); %get parameters
n = par.n; % number of segments

% relative angles
phi = sym('phi', [n,1]);
phi_t = phi;
for i = 1:n
    phi_t(i) = str2sym(['phi' num2str(i) '(t)']);
end
% relative angular velocities and accelerations
phid = sym('phid', [n,1]);
phid_t = diff(phi_t);
phidd = sym('phidd', [n,1]);
phidd_t = diff(phid_t);
% center of mass coordinates
syms L
x(1) = L/2*sin(phi_t(1));
for i = 2:n
    x(i) = 0;
    for j = 1:i-1
        x(i) = x(i) + L*sin(sum(phi_t(1:j)));
    end
    x(i) = x(i) + L/2*sin(sum(phi_t(1:i)));
end
x = x.';
y(1) = -L/2*cos(phi_t(1));
for i = 2:n
    y(i) = 0;
    for j = 1:i-1
        y(i) = y(i) - L*cos(sum(phi_t(1:j)));
    end
    y(i) = y(i) - L/2*cos(sum(phi_t(1:i)));
end
y = y.';
% velocities
xd = diff(x);
yd = diff(y);

% kinetic energy
syms m theta
for i = 1:n
    omega(i) = sum(phid_t(1:i));
end
omega = omega.';
T = 0;
for i = 1:n
    T = T + 1/2*m*(xd(i)^2 + yd(i)^2) + 1/2*theta*omega(i)^2;
end

%potential energy
kt = sym('kt', [1, n+1]);
U = 0;
for i = 1:n
    U = U + 1/2*kt(i)*(phi_t(i)-par.phi_r(i))^2;
end
U = U + 1/2*kt(end)*(pi-sum(phi_t)-par.phi_r(end))^2;

% dissipative potential
b = sym('b', [1, n+1]);
D = 0;
for i = 1:n
    D = D + 1/2*b(i)*phid_t(i)^2;
end
D = D + 1/2*b(end)*(-sum(phid_t))^2;

% derivatives
Tsub = subs(T, phid_t, phid);
Tphid = jacobian(Tsub, phid).';
Tphid = subs(Tphid, phid, phid_t);
Tphid = diff(Tphid);
Tphid = subs(Tphid, [phi_t; phid_t; phidd_t], [phi; phid; phidd]);

Tsub = subs(T, [phi_t; phid_t], [phi; phid]);
Tphi = jacobian(Tsub, phi).';

Dsub = subs(D, phid_t, phid);
Dphid = jacobian(Dsub, phid).';

Usub = subs(U, phi_t, phi);
Uphi = jacobian(Usub, phi).';

% Lagrange equation
Lg = Tphid - Tphi + Dphid + Uphi;
% sort out the coeffitients
for i = 1:n
    coef(i,:) = coeffs(Lg(i), phidd);
end
coef = flip(coef,2);
M = coef(1:n,1:n);
K = coef(:,end);

% constraints
%  http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node90.html
syms Diam
xend = 0;
yend = 0;
for i = 1:n
    xend = xend + L*sin(sum(phi(1:i)));
    yend = yend - L*cos(sum(phi(1:i)));
end
c = [xend - Diam; yend];   %h(q) = 0

C = jacobian(c,phi);
H = jacobian(C*phid, phi)*phid;

% contact force calculation
syms h d epsilon v mu
for k = 1:n   % external forces and positions where they're applied
    % position of contact point
    Xc(k) = x(k) + h*sin(sum(phi(1:k)) - pi/2);
    Yc(k) = y(k) - h*cos(sum(phi(1:k)) - pi/2);
    delta(k) = Yc(k) + Diam/2 + h - d; % distance from contact surface
    Fy(k) = epsilon * exp(-delta(k)/epsilon); % global normal contact force
    Fx(k) = -mu*Fy(k)*tanh(10*(diff(Xc(k))-v)); % global friction force
end
% get rid of (t)
Xc = subs(Xc.', [phi_t; phid_t], [phi; phid]);
Yc = subs(Yc.', [phi_t; phid_t], [phi; phid]);
delta = subs(delta.', [phi_t; phid_t], [phi; phid]);
Fx = subs(Fx.', [phi_t; phid_t], [phi; phid]);
Fy = subs(Fy.', [phi_t; phid_t], [phi; phid]);

Q = sym('Q', [n,1]);
for j = 1:n % generalized forces
    Q(j) = 0;
    for k = 1:n
        Q(j) = Q(j) + [Fx(k), Fy(k)] * jacobian([Xc(k),Yc(k)],phi(j));
    end
end

toc
%% substitute constants
tic
% substitute the numerical values
D_num = par.D;
m_num = par.m;
L_num = par.L;
k_num = par.k;
b_num = par.b;
theta_num = par.theta;
h_num = par.h;
mu_num = par.mu;
d_num = par.d;

epsilon_num = par.epsilon;

params = [Diam;m;L;kt.';b.';theta;h;mu;epsilon];
params_num = [D_num;m_num;L_num;k_num';b_num';theta_num;...
    h_num;mu_num;epsilon_num];
vd_sym = [v;d];

Msub = subs(M,params,params_num);
Csub = subs(C,params,params_num);
Ksub = subs(K,params,params_num);
Qsub = subs(Q,params,params_num);
Hsub = subs(H,params,params_num);
Fxsub = subs(Fx,params,params_num);
Fysub = subs(Fy,params,params_num);

Mfun = matlabFunction(Msub);
Cfun = matlabFunction(Csub);
Kfun = matlabFunction(Ksub);
Qfun = matlabFunction(Qsub);
Hfun = matlabFunction(Hsub);
Fxfun = matlabFunction(Fxsub);
Fyfun = matlabFunction(Fysub);

save('eq_of_motion_data.mat','Mfun','Cfun','Kfun','Qfun','Hfun','Fxfun');
toc
%% solve ODE
tic
% calculate static deformed shape
% use relaxed state as initial condition
x0 = zeros(1,(n+1)*3);
x0(1:n+1) = par.phi_r;
options = optimoptions('fsolve','MaxFunctionEvaluations',80000,...
    'MaxIterations', 5000);
[x, fval, exitflag] = fsolve(@static_equations_approx, x0, options);
if exitflag < 0
    error('No initial deformed shape found!');
end
phi0 = x(1:n);
toc

% modify IC
% phi0(1) = phi0(1)+0.05;
% init = get_dynamic_IC(phi0);
load init.mat
if par.force_mode == 1
    init = [init, 0, 0];
end
animateSensor(0,[init(1:n),x(n+1:end)]); title('initial shape');
% dynamic simulation
tspan = 0:0.01:5;
options = odeset('RelTol',1e-6,'AbsTol',1e-8, 'BDF', 'on');
[t,Y] = ode15s(@eq_of_motion,tspan,init');
toc
%% animation
% Calculating joint coordinates for animation purposes
x = L_num*sin(Y(:,1));
y = -L_num*cos(Y(:,1));
if n > 1
    for i = 2:n
        x(:,i) = x(:,i-1) + L_num*sin(sum(Y(:,1:i),2));
        y(:,i) = y(:,i-1) - L_num*cos(sum(Y(:,1:i),2));
    end
end

ang = Y(:,1:n);

% Set up first frame
figure('Color', 'white')%,'units','normalized','position',[0 0 1 1])

subplot(2,1,1)
cmap = colormap(jet(n));
for i = 1:n
    plot(t, ang(:,i), 'LineWidth', 2, 'Color', cmap(i,:))
    hold on
end

for i = 1:n
    hh1(i) = line(t(1), ang(1,i), 'Marker', '.', 'MarkerSize', 20, ...
        'Color', 'r');
end
xlabel('t [s]')
ylabel('relative angle [rad]')

subplot(2,1,2)
xplot = [0, x(1,1)];
yplot = [0, y(1,1)];
if n > 1
    for i = 2:n
        xplot(i,:) = [xplot(i-1,2), x(1,i)];
        yplot(i,:) = [yplot(i-1,2), y(1,i)];
    end
end
% contact surface
hh3 = line([-0.5*D_num 1.5*D_num],...
    [-D_num/2-h_num+d_num(t(1)) -D_num/2-h_num+d_num(t(1))],'LineWidth',1,'Color','k');
hold on
% top of sensor
line([0,D_num], [0,0],'LineWidth',1,'Color','k');
for i = 1:n
    hh2(i) = plot(xplot(i,:), yplot(i,:), '.-', 'MarkerSize', 20, 'LineWidth', 2, ...
        'Color', cmap(i,:));
end

hold off
axis equal
axis([-0.5*D_num 1.5*D_num -D_num 0.5*D_num])
ht = title(sprintf('t = %0.2f s', t(1)));


% Get figure size
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);


% % Preallocate data (for storing frame data)
% mov = zeros(height, width, 1, length(t), 'uint8');
%
% Loop through by changing XData and yData
for id = 1:length(t)
    % Update graphics data. this is more efficient than recreating plots.
    for j = 1:n
        set(hh1(j), 'XData', t(id), 'yData', ang(id, j))
        xplot = [0, x(id,1)];
        yplot = [0, y(id,1)];
        set(hh2(1), 'XData', xplot, 'yData', yplot)
        for i = 2:n
            xplot(i,:) = [xplot(i-1,2), x(id,i)];
            yplot(i,:) = [yplot(i-1,2), y(id,i)];
            set(hh2(i), 'XData', xplot(i,:), 'yData', yplot(i,:))
        end
        set(ht, 'String', sprintf('t = %0.2f s', t(id)))
        set(hh3, 'XData', [-0.5*D_num 1.5*D_num],...
            'YData', [-D_num/2-h_num+d_num(t(id)) -D_num/2-h_num+d_num(t(id))]);
    end
    drawnow;
    pause(0.03)
end
%
%     % Get frame as an image
%     f = getframe(gcf);
%
%     % Create a colormap for the first frame. For the rest of the frames,
%     % use the same colormap
%     if id == 1
%         [mov(:,:,1,id), map] = rgb2ind(f.cdata, 256, 'nodither');
%     else
%         mov(:,:,1,id) = rgb2ind(f.cdata, map, 'nodither');
%     end
% end

%
% % Create animated GIF
% imwrite(mov, map, 'animation.gif', 'Delaytime', 0, 'LoopCount', inf)

%% calculate forces and friction coeffitient
Fxval = zeros(length(t),par.n);
Fyval = zeros(length(t),par.n);
for i = 1:length(t)
    switch par.force_mode
        case 0
            Fxarg = num2cell([par.d(t(i)), Y(i,:), par.v(t(i))]);
            Fyarg = num2cell([par.d(t(i)), Y(i,1:par.n)]);
        case 1
            Fxarg = num2cell([par.d(t(i)), Y(i,1:2*par.n), Y(i,end)]);
            Fyarg = num2cell([par.d(t(i)), Y(i,1:par.n)]);
    end
    
    Fxval(i,:) = Fxfun(Fxarg{:})';
    Fyval(i,:) = Fyfun(Fyarg{:})';
    
    Fxsum(i) = sum(Fxval(i,:));
    Fysum(i) = sum(Fyval(i,:));
end
figure
plot(t,Fxsum, 'LineWidth', 2);
xlabel('t [s]'); ylabel('F_f [10^{-3} N]', 'Interpreter',  'tex');
figure
plot(t,Fysum, 'LineWidth', 2);
xlabel('t [s]'); ylabel('F_n [10^{-3} N]', 'Interpreter',  'tex');

%% plots to save

% relative angles
figure
for i = 1:n
    plot(t, ang(:,i), 'LineWidth', 2, 'Color', cmap(i,:))
    hold on
end
xlabel('t [s]')
ylabel('relative angle [rad]')

% relative angular velocities
figure
angvel = Y(:,n+1:2*n);
for i = 1:n
    plot(t, angvel(:,i), 'LineWidth', 2, 'Color', cmap(i,:))
    hold on
end
xlabel('t [s]')
ylabel('relative anglular velocity [rad/s]')

% contact depth
figure
plot(t, par.d(t),'LineWidth', 2)
xlabel('t [s]')
ylabel('d [mm]')

% relative velocity
figure
switch par.force_mode
    case 0
        plot(t, par.v(t),'LineWidth', 2)
    case 1
        plot(t, Y(:,end),'LineWidth', 2)
end
xlabel('t [s]')
ylabel('v [mm/s]')

%% sensor shape
figure
k = find(t == 2.5);
xplot = [0, x(k,1)];
yplot = [0, y(k,1)];
if n > 1
    for i = 2:n
        xplot(i,:) = [xplot(i-1,2), x(k,i)];
        yplot(i,:) = [yplot(i-1,2), y(k,i)];
    end
end
% contact surface
hh3 = line([-0.5*D_num 1.5*D_num],...
    [-D_num/2-h_num+d_num(t(k)) -D_num/2-h_num+d_num(t(k))],'LineWidth',1,'Color','k');
hold on
% top of sensor
line([0,D_num], [0,0],'LineWidth',1,'Color','k');
% segments
for i = 1:n
    hh2(i) = plot(xplot(i,:), yplot(i,:), '.-', 'MarkerSize', 20, 'LineWidth', 2, ...
        'Color', cmap(i,:));
end

hold off
axis equal
axis([-0.5*D_num 1.5*D_num -D_num 0.5*D_num])
ht = title(sprintf('t = %0.2f s', t(k)),'FontSize',16);
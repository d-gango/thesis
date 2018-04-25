%% 
clear all
% dof
n = 3;
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
% velovities
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
syms g
U = 0;
for i = 1:n
    U = U + m*g*y(i);
end

% dissipative potential
syms b
D = 0;
for i = 1:n
    D = D + 1/2*b*omega(i)^2;
end

% derivatives
Tsub = subs(T, phid_t, phid);
Tphid = jacobian(Tsub, phid).';
Tphid = subs(Tphid, phid, phid_t);
Tphid = diff(Tphid);
Tphid = simplify(subs(Tphid, [phi_t; phid_t; phidd_t], [phi; phid; phidd]));

Tsub = subs(T, [phi_t; phid_t], [phi; phid]);
Tphi = simplify(jacobian(Tsub, phi).');

Dsub = subs(D, phid_t, phid);
Dphid = simplify(jacobian(Dsub, phid).');

Usub = subs(U, phi_t, phi);
Uphi = simplify(jacobian(Usub, phi).');

% Lagrange equation
Lg = Tphid - Tphi + Dphid + Uphi;
% sort out the coeffitients
for i = 1:n
    coef(i,:) = coeffs(Lg(i), phidd);
end
coef = flip(coef,2);
M = coef(1:n,1:n);
C = coef(:,end);

%%
% substitute the numerical values
m_num = 1;
L_num = 1;
g_num = 9.81;
b_num = 0;
theta_num = 1/12*m_num*L_num^2;

M = subs(M, [m,L,g,b,theta], [m_num,L_num,g_num,b_num,theta_num]);
C = subs(C, [m,L,g,b,theta], [m_num,L_num,g_num,b_num,theta_num]);

phid_t = phi;
for i = 1:n
    phid_t(i) = str2sym(['phid' num2str(i) '(t)']);
end
M = subs(M, [phi; phid], [phi_t; phid_t]);
C = subs(C, [phi; phid], [phi_t; phid_t]);
acc = M \ (-C);

eq = [phid_t; acc];
f = odeFunction(eq, [phi_t;phid_t]);

%% solve ODE
init = [pi/3, pi/4 -pi/2, 0, 0, 0]';
tspan = linspace(0, 30, 600);

[t,Y] = ode45(f,tspan,init);

%% animation
L = 1;
% Calculating joint coordinates for animation purposes
x = L*sin(Y(:,1));
y = -L*cos(Y(:,1));
if n > 1
    for i = 2:n
        x(:,i) = x(:,i-1) + L*sin(sum(Y(:,1:i),2));
        y(:,i) = y(:,i-1) - L*cos(sum(Y(:,1:i),2));
    end
end

ang = Y(:,1:n);

% Set up first frame
figure('Color', 'white')
subplot(2,1,1)
plot(t, ang, 'LineWidth', 2)
colors = 'brygkmc';
for i = 1:n
    hh1(i) = line(t(1), ang(1,i), 'Marker', '.', 'MarkerSize', 20, ...
        'Color', colors(i));
end
xlabel('time (sec)')
ylabel('angle (rad)')

subplot(2,1,2)
xplot = [0, x(1,1)];
yplot = [0, y(1,1)];
for i = 2:n
    xplot(i,:) = [xplot(i-1,2), x(1,i)];
    yplot(i,:) = [yplot(i-1,2), y(1,i)];
end
for i = 1:n
    hh2(i) = plot(xplot(i,:), yplot(i,:), '.-', 'MarkerSize', 20, 'LineWidth', 2);
    hold on
end
hold off
axis equal
axis([-n*L n*L -n*L n*L])
ht = title(sprintf('time: %0.2f sec', t(1)));

% Get figure size
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);

% % Preallocate data (for storing frame data)
% mov = zeros(height, width, 1, length(t), 'uint8');

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
        set(ht, 'String', sprintf('time: %0.2f sec', t(id)))
    end
    drawnow;
    pause(0.03)
end

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
%% check total energy
phid_t = diff(phi_t);
syms L
kin = subs(T, [m,L,g,theta], [m_num,L_num,g_num,theta_num]);
pot = subs(U, [m,L,g,theta], [m_num,L_num,g_num,theta_num]);
total_energy = zeros(length(t),1);
for i = 1:length(t)
    total_energy(i) = double(subs(kin+pot, [phi_t; phid_t], Y(i,:)'));
end
energy_fluctuation = max(total_energy) - min(total_energy);
rel_energy_fluctuation = energy_fluctuation / mean(total_energy);
disp(['Energy fluctuation: ', num2str(rel_energy_fluctuation*100), ' %'])

figure
plot(t, total_energy)
title('Total energy')
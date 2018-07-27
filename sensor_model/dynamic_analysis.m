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
    delta(k) = Yc(k) + Diam/2 - d; % distance from contact surface
    Fy(k) = epsilon * exp(-delta(k)/epsilon); % global normal contact force
    Fx(k) = mu*Fy(k)*tanh(10*(v-diff(Xc(k)))); % global friction force
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
%     for k = 1:n
%         Q(j) = Q(j) + [Fx(k), Fy(k)] * jacobian([Xc(k),Yc(k)],phi(j));
%     end
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

Mfun = matlabFunction(Msub);
Cfun = matlabFunction(Csub);
Kfun = matlabFunction(Ksub);
Qfun = matlabFunction(Qsub);
Hfun = matlabFunction(Hsub);

save('eq_of_motion_data.mat','Mfun','Cfun','Kfun','Qfun','Hfun');
toc
%% system linearisation
% stepsize
h = 0.000001;
% equilibrium
x0 = [par.phi_r(1:n), zeros(1,n)];

Alin = zeros(2*n);
% calculate partial derivatives
for i = 1:2*n
    x = x0.';
    x(i) = x(i)+h;
    x0plush = eq_of_motion_without_forces(0,x);
    x = x0.';
    x(i) = x(i)-h;
    x0minush = eq_of_motion_without_forces(0,x);
    dfdx = (x0plush - x0minush) ./ (2*h);
    Alin(:,i) = dfdx;
end

% eigenvalues and vectors
[V,D] = eig(Alin);

% modes
modes = inv(V);

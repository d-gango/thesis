%% equation of motion
tic
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
syms k
U = 0;
for i = 1:n
    if i == 1
        U = U + k*phi_t(i)^2;  % double spring stiffness
    else
    U = U + 1/2*k*phi_t(i)^2;
    end
end
U = U + k*(pi-sum(phi_t))^2;  % double spring stiffness

% dissipative potential
syms b
D = 0;
for i = 1:n
    D = D + 1/2*b*phid_t(i)^2;
end
D = D + 1/2*b*(-sum(phid_t))^2;

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

% constraints 
%  http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node90.html
syms Diam
xend = 0;
yend = 0;
for i = 1:n
    xend = xend + L*sin(sum(phi(1:i)));
    yend = yend - L*cos(sum(phi(1:i)));
end
h = [xend - Diam; yend];   %h(q) = 0
Cm = jacobian(h, [phi;phid]);
toc
%% substitute constants and generate function handle
% substitute the numerical values
D = 40;
m_num = 100/n;
L_num = D*sin(pi/(2*n));
k_num = 10000;
b_num = 10;
theta_num = 1/12*m_num*L_num^2;

M = subs(M, [m,L,k,b,theta], [m_num,L_num,k_num,b_num,theta_num]);
C = subs(C, [m,L,k,b,theta], [m_num,L_num,k_num,b_num,theta_num]);
Cm = subs(Cm, [L, Diam], [L_num, D]);

phid_t = phi;
for i = 1:n
    phid_t(i) = str2sym(['phid' num2str(i) '(t)']);
end
% M = subs(M, [phi; phid], [phi_t; phid_t]);
% C = subs(C, [phi; phid], [phi_t; phid_t]);
% Cm = subs(Cm, [phi; phid], [phi_t; phid_t]);
acc = M \ (-C);
fq = [phid_t; acc];
eq = fq + (Cm.')/(Cm*Cm.') * (-Cm*fq);
% linearization
Alin = jacobian(eq,[phi; phid]);
% eqiulibrium state
phi0 = ones(n,1)*pi/n;
phi0(1) = phi0(1)/2;
phid0 = zeros(n,1);
% substitute equilibrium state
Alin = double(subs(Alin, [phi;phid], [phi0;phid0]));
% eigenvalues
[V,D] = eig(Alin);
lambda = diag(D)
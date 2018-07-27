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
h = [xend - Diam; yend];   %h(q) = 0

C = jacobian(h,phi);
H = jacobian(C*phid, phi)*phid;

% first order ODE eq of motion
Mbar = sym('m',[2*n,2*n]);
Mbar(1:n,n+1:2*n) = zeros(n,n);
Mbar(n+1:2*n,1:n) = zeros(n,n);
Mbar(1:n,1:n) = eye(n);
Mbar(n+1:2*n,n+1:2*n) = M;
kbar = [-phid;K];
Q = zeros(n,1);
Qbar = zeros(2*n,1);
Cbar = [zeros(n,length(h)); C.'];



% constraint enforcement !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
f = C*phid;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
A1 = C*(M\(C.'));
A2 = (M\(-K+Q));
A3 = (C*A2);
A4 = -H-f-A3;

% for simlicity, phid = 0 when liearising around a fixed point
% A1 = subs(A1,phid, zeros(size(phid)));
% A4 = subs(A4,phid, zeros(size(phid)));

lambda = A1\A4;

toc
%% substitute constants and generate function handle
%the numerical values
D_num = 40;
m_num = 100/n;
L_num = D_num*sin(pi/(2*n));
k_num = 10000;
b_num = 1000;
theta_num = 1/12*m_num*L_num^2;

% linearized system matrix
dkdq = jacobian(K,phi);
dkdqd = jacobian(K,phid);
dQdq = jacobian(Q,phi);
dcldq = jacobian(C.' * lambda,phi);
dcldqd = jacobian(C.' * lambda,phid);
% eqiulibrium state
phi0 = ones(n,1)*pi/n;
phi0(1) = phi0(1)/2;
phid0 = zeros(n,1);
% substitute equilibrium state and numerical values
M_num = double( subs(M, [Diam,m,L,k,b,theta,phi.',phid.'],...
               [D_num,m_num,L_num,k_num,b_num,theta_num,phi0.',phid0.']));
dkdq_num = double( subs(dkdq, [Diam,m,L,k,b,theta,phi.',phid.'],...
               [D_num,m_num,L_num,k_num,b_num,theta_num,phi0.',phid0.']));
dkdqd_num = double( subs(dkdqd, [Diam,m,L,k,b,theta,phi.',phid.'],...
               [D_num,m_num,L_num,k_num,b_num,theta_num,phi0.',phid0.']));
dQdq_num = double( subs(dQdq, [Diam,m,L,k,b,theta,phi.',phid.'],...
               [D_num,m_num,L_num,k_num,b_num,theta_num,phi0.',phid0.']));
dcldq_num = double( subs(dcldq, [Diam,m,L,k,b,theta,phi.',phid.'],...
               [D_num,m_num,L_num,k_num,b_num,theta_num,phi0.',phid0.']));
dcldqd_num = double( subs(dcldqd, [Diam,m,L,k,b,theta,phi.',phid.'],...
               [D_num,m_num,L_num,k_num,b_num,theta_num,phi0.',phid0.']));

A = [zeros(n,n), eye(n);...
    M_num\(-dkdq_num + dQdq_num + dcldq_num) , M_num\(-dkdqd_num + dcldqd_num)];
% eigenvalues
[V,D] = eig(A);
eigenvalues = diag(D)
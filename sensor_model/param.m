function par = param()
par.n = 20;  % number of segments
par.epsilon = 0.06; % parameter for contact force approximation
par.v = @(t) 0;  % relative velocity of contact surface
par.d = @(t) 6;  % contact depth

par.batch = 0; % set to 1 for batch run
par.offset = 0; % springs relaxed in equilibrium

par.D = 40;  % sensor diameter
par.m = 100/par.n;  % mass of one segment
par.L = par.D*sin(pi/(2*par.n)); % length of one segment
par.k = ones(1,par.n+1)*10000; % spring stiffness
par.k(1) = 2*par.k(1); par.k(end) = par.k(1);
par.b = ones(1,par.n+1)*1000; % damping coeff.
par.b(1) = 2*par.b(1); par.b(end) = par.b(1);
par.theta = 1/12*par.m*par.L^2; % moment of inertia
par.h = 0.3;  % length of perpendicular contact part
par.mu = 0.7; % sliding friction coeff.

phi_r = ones(1,par.n+1)*pi/par.n;  % relaxed state angles
phi_r(1) = phi_r(1)/2;
phi_r(end) = phi_r(end)/2;
par.phi_r = phi_r;

% new kt calculation
par.t = 0.3; % thickness [mm]
par.E = 25000; % Young's modulus [kPa]
psi_r = getPsi(par.phi_r);
dy = par.L*cos(psi_r);
kt =(par.E * par.t^3 / 6).*...
    1./(atan(dy(1:par.n/2) ./ sqrt((par.D/2)^2 - dy(1:par.n/2).^2)));
par.k = [kt, 10e6, flip(kt)];
end
%===================================================================
function psi = getPsi(phi)
n = length(phi);
psi = zeros(1,n);
for j = 1:n
    psi(j) = sum(phi(1:j));
end
end
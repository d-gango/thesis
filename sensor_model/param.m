function par = param()
par.n = 20;  % number of segments
par.epsilon = 0.06; % parameter for contact force approximation
par.v = @(t) 0;  % relative velocity of contact surface
par.d = @(t) 5;  % contact depth

par.batch = 0; % set to 1 for batch run
par.offset = 1; % springs relaxed in equilibrium

par.D = 39.5;  % sensor diameter
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

% % new kt calculation
% par.t = 1; % thickness [mm]
% par.E = 25000; % Young's modulus [kPa]
% psi_r = getPsi(par.phi_r);
% y = zeros(1,par.n+1);
% for i = 2:par.n+1
%     y(i) = y(i-1) - par.L*cos(psi_r(i-1));
% end
% d = zeros(1,par.n/2);
% for i=1:par.n/2
%     d(i) = 2*pi*sqrt(par.D^2/4 - ((y(i)+y(i+1))/2)^2);
% end
% kt = (par.E .* d .* par.t^3) ./ (12*par.L);
% par.k = [kt, kt(end), flip(kt)];
end
%===================================================================

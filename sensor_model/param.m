function par = param()
par.n = 10;  % number of segments
par.epsilon = 0.06; % parameter for contact force approximation
par.v = @(t) heaviside(t-1);  % relative velocity of contact surface
par.d = @(t) 5;  % contact depth

par.batch = 0; % set to 1 for batch run
par.offset = 1; % springs relaxed in equilibrium

par.p =0; % pressure inside
par.D = 38.5;  % sensor diameter
par.m = 100/par.n;  % mass of one segment
par.L = par.D*sin(pi/(2*par.n)); % length of one segment
par.k = ones(1,par.n+1)*50000; % spring stiffness
par.k(1) = 2*par.k(1); par.k(end) = par.k(1);
par.b = ones(1,par.n+1)*1000; % damping coeff.
par.b(1) = 2*par.b(1); par.b(end) = par.b(1);
par.theta = 1/12*par.m*par.L^2; % moment of inertia
par.h = 0.5;  % length of perpendicular contact part
par.mu = 0.5; % sliding friction coeff.
par.c = 0;%ones(1,par.n) * 1/10000;

phi_r = ones(1,par.n+1)*pi/par.n;  % relaxed state angles
phi_r(1) = phi_r(1)/2;
phi_r(end) = phi_r(end)/2;
par.phi_r = phi_r;
psi_r = getPsi(par.phi_r);
% surface of segments
par.A = par.D/2 * par.L * pi * abs(cos(psi_r(1:par.n)));

% % new kt calculation
% par.t = 1; % thickness [mm]
% par.E = 45000; % Young's modulus [kPa]
% psi_r = getPsi(par.phi_r);
% y = zeros(1,par.n+1);
% for i = 2:par.n+1
%     y(i) = y(i-1) - par.L*cos(psi_r(i-1));
% end
% % two different methods for width, CHOOSE ONE!
% d = zeros(1,par.n/2);  % simple mean
% w = zeros(1,par.n/2);  % integral mean
% for i=1:par.n/2
%     w(i) =2*pi/par.L * integral(@(ybar) sqrt(par.D^2/4 - (y(i) + ybar.*cos(psi_r(i))).^2),0,par.L);
%     d(i) = 2*pi*sqrt(par.D^2/4 - ((y(i)+y(i+1))/2)^2);
% end
% kt = (par.E .* w .* par.t^3) ./ (12*par.L);
% par.k = [kt, kt(end), flip(kt)];
% % c = par.L ./ (w .* par.t .* par.E);
% % par.c = [c, flip(c)];
% end
%===================================================================

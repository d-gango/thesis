function par = param()
par.n = 20;  % number of segments
par.epsilon = 0.06; % parameter for contact force approximation
par.v = @(t) 10;  % relative velocity of contact surface
par.d = @(t) 6;  % contact depth

par.batch = 1; % set to 1 for batch run

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
end
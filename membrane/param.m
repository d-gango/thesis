function par = param()
par.n = 4;
par.epsilon = 0.05;
par.v = @(t) 0;
par.d = @(t) 5;

par.D = 40;
par.m = 100/par.n;
par.L = par.D*sin(pi/(2*par.n));
par.k = 10000;
par.b = 1000;
par.theta = 1/12*par.m*par.L^2;
par.h = 0.3;
par.mu = 1;

phi_r = ones(1,par.n+1)*pi/par.n;
phi_r(1) = phi_r(1)/2;
phi_r(end) = phi_r(end)/2;
par.phi_r = phi_r;
end
function xdot = odefun(t,x)
% load symbolic matrices
load('eq_of_motion_data.mat');
% substitute numeric values
C = double(subs(C, [params; phi; phid], [params_num; x]));
Cbar = double(subs(Cbar, [params; phi; phid], [params_num; x]));
K = double(subs(K, [params; phi; phid], [params_num; x]));
kbar = double(subs(kbar, [params; phi; phid], [params_num; x]));
M = double(subs(M, [params; phi; phid], [params_num; x]));
Mbar = double(subs(Mbar, [params; phi; phid], [params_num; x]));
Q = double(subs(Q, [params; phi; phid], [params_num; x]));
Qbar = double(subs(Qbar, [params; phi; phid], [params_num; x]));
H = double(subs(H, [params; phi; phid], [params_num; x]));
% calculate xdot
lambda = (C*(M\(C.')))\(-H-(C*(M\(-K+Q))));
xdot = Mbar\(-kbar + Qbar + Cbar*lambda);
end
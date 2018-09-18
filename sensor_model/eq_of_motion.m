function xdot = eq_of_motion(t,x)
% load symbolic matrices
load('eq_of_motion_data.mat');
par = param();
phicell = num2cell(x(1:par.n));
xcell = num2cell(x(1:2*par.n));
Qcell = num2cell([par.d(t), x(1:par.n).',x(2*par.n+1:end).' , par.v(t)]);
phi_bardotcell = num2cell([x(2*par.n+1:end).', par.v(t)]);
% substitute numeric values
C = Cfun(phicell{:});
K = Kfun(xcell{:});
M = Mfun(phicell{:});
Q = Qfun(Qcell{:});
H = Hfun(xcell{:});
phi_bardot = phi_bardotfun(phi_bardotcell{:});
% calculate xdot
lambda = (C*(M\(C.')))\(-H-(C*(M\(-K+Q))));
xdot = [x(par.n+1:2*par.n); M\(-K + Q + C.'*lambda); phi_bardot];
disp(t)
end
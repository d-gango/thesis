function xdot = eq_of_motion_without_forces(t,x)
% load symbolic matrices
load('eq_of_motion_data.mat');
par = param();
phicell = num2cell(x(1:par.n));
xcell = num2cell(x);
% substitute numeric values
C = Cfun(phicell{:});
K = Kfun(xcell{:});
M = Mfun(phicell{:});
Q = Qfun();
H = Hfun(xcell{:});
% calculate xdot
lambda = (C*(M\(C.')))\(-H-(C*(M\(-K+Q))));
xdot = [x(par.n+1:end); M\(-K + Q + C.'*lambda)];
end
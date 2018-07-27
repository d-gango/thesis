function xdot = eq_of_motion(t,x)
% load symbolic matrices
load('eq_of_motion_data.mat');
par = param();
phicell = num2cell(x(1:par.n));
xcell = num2cell(x);
dxvcell = num2cell([par.d(t), x.', par.v(t)]);
% substitute numeric values
C = Cfun(phicell{:});
K = Kfun(xcell{:});
M = Mfun(phicell{:});
Q = Qfun(dxvcell{:});
H = Hfun(xcell{:});
% calculate xdot
lambda = (C*(M\(C.')))\(-H-(C*(M\(-K+Q))));
xdot = [x(par.n+1:end); M\(-K + Q + C.'*lambda)];
disp(t)
end
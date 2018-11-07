function xdot = eq_of_motion(t,x)
global Cfun Kfun Mfun Qfun Hfun
par = param();
phicell = num2cell(x(1:par.n));
xcell = num2cell(x(1:2*par.n));
switch par.force_mode % force_mode -> we simulate the surface as a moving body
    case 0
        dxvcell = num2cell([par.d(t), x(1:2*par.n).', par.v(t)]);
    case 1
        dxvcell = num2cell([par.d(t), x(1:2*par.n).', x(end-1)]);
        fxcell = num2cell([par.d(t), x(1:2*par.n).', x(end-1)]);
end
% substitute numeric values
C = Cfun(phicell{:});
K = Kfun(xcell{:});
M = Mfun(phicell{:});
Q = Qfun(dxvcell{:});
H = Hfun(xcell{:});
if par.force_mode == 1
    Fx = Fxfun(fxcell{:});
    Fxsum = sum(Fx);
end

% calculate xdot
lambda = (C*(M\(C.')))\(-H-(C*(M\(-K+Q))));
switch par.force_mode
    case 0
        xdot = [x(par.n+1:end); M\(-K + Q + C.'*lambda)];
    case 1
        xdot = [x(par.n+1:2*par.n); M\(-K + Q + C.'*lambda); x(end); (par.F_ext-Fxsum)/par.m_surf];
end

%disp(t)
end
function xdot = eq_of_motion_spinodal(t,x)
% x = [phi, phidot, phi_bar]
% force mode x = x = [phi, phidot, phi_bar ,x, xdot]
global Cfun Kfun Mfun Qfun Hfun phi_bardotfun Fxfun
par = param();
switch par.force_mode
    case 0
        phicell = num2cell(x(1:par.n));
        xcell = num2cell(x(1:2*par.n));
        Qcell = num2cell([par.d(t), x.' , par.v(t)]);
        phi_bardotcell = num2cell([x.', par.v(t)]);
    case 1
        phicell = num2cell(x(1:par.n));
        xcell = num2cell(x(1:2*par.n));
        Qcell = num2cell([par.d(t), x(1:3*par.n).' , x(end)]);
        phi_bardotcell = num2cell([x(1:3*par.n).', x(end)]);
        Fxcell = num2cell([par.d(t), x(1:3*par.n).', x(end)]);
end

% substitute numeric values
C = Cfun(phicell{:});
K = Kfun(xcell{:});
M = Mfun(phicell{:});
Q = Qfun(Qcell{:});
H = Hfun(xcell{:});
phi_bardot = phi_bardotfun(phi_bardotcell{:});
if par.force_mode == 1
    Fx = Fxfun(Fxcell{:});
    Fxsum = sum(Fx); % net friction force
end

% calculate xdot
lambda = (C*(M\(C.')))\(-H-(C*(M\(-K+Q))));
switch par.force_mode
    case 0
    xdot = [x(par.n+1:2*par.n); M\(-K + Q + C.'*lambda); phi_bardot];
    case 1
    % net horizontal force acting on the surface
    Fsurf = par.F_ext(t) - Fxsum;
    xdot = [x(par.n+1:2*par.n); M\(-K + Q + C.'*lambda); phi_bardot; x(end); Fsurf/par.m_surf];
end
    
%disp(t)
end
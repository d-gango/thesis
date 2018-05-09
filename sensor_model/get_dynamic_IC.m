% construct initial conditions consistent with constraints
function init = get_dynamic_IC(phi0)
par = param();
n = par.n;
phi = sym('phi', [n,1]);
syms D L
par = param();
D_num = par.D;
L_num = par.L;
xend = 0;
yend = 0;
for i = 1:n
    xend = xend + L*sin(sum(phi(1:i)));
    yend = yend - L*cos(sum(phi(1:i)));
end
eqs = [xend - D; yend];
eqs = subs(eqs, [D,L,phi(1:end-2).'], [D_num,L_num,phi0(1:end-2)]);
sol = vpasolve(eqs, phi(end-1:end), phi0(end-1:end));
if isempty(struct2array(sol))
    error('No initial condition found.')
end
phi_init = double([phi0(1:end-2), struct2array(sol)]);
phid_init = zeros(1,n);

init = [phi_init, phid_init];
end
function init = findIC(n,phi0)
phi = sym('phi', [n,1]);
syms D L
D_num = 40;
L_num = D_num*sin(pi/(2*n));
xend = 0;
yend = 0;
for i = 1:n
    xend = xend + L*sin(sum(phi(1:i)));
    yend = yend - L*cos(sum(phi(1:i)));
end
eqs = [xend - D; yend];
eqs = subs(eqs, [D,L,phi(1:end-2).'], [D_num,L_num,phi0]);
sol = vpasolve(eqs, phi(end-1:end), [pi/n pi/n]);
if isempty(struct2array(sol))
    error('No initial condition found.')
end
phi_init = double([phi0, struct2array(sol)]);
phid_init = zeros(1,n);

init = [phi_init, phid_init];
end
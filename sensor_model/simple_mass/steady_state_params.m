a = 0.0349;
bb = 0.0489;
c = 1e-3;
R = 100;
V_star = 1e-6;
mu_star = 0.369;

phi_ss = @(vbar)(1+R)./(1+R.*vbar);
mufun = @(vbar) a.*asinh(vbar./2 .* exp((mu_star + bb.*log(c+(1+R)./(1+R.*vbar)))./a));

v = logspace(-5, 5,1000);
mu = mufun(v);
phi = phi_ss(v);

figure
semilogx(v,mu,'LineWidth',2)
xlabel('$\frac{V}{V_*}$','Interpreter','latex','FontSize',18);
ylabel('$\mu_{ss}$','Interpreter','latex','FontSize',18);

figure
semilogx(v,phi,'LineWidth',2)
xlabel('$\frac{V}{V_*}$','Interpreter','latex','FontSize',18);
ylabel('$\phi_{ss}$','Interpreter','latex','FontSize',18);

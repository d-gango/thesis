% !!! rewritten for dynamic simulation !!!!
function eq = static_equations_approx(x)
par = param();
% numeric parameters
n = par.n;
D = par.D;
if par.batch
    d = getGlobald();
else
    d = par.d(0);
end
L = par.L;
kt = par.k;
h = par.h;
mu = par.mu;
phi_r = par.phi_r;
epsilon = par.epsilon;
v = par.v;
eq = [];
% relative angles
phivector = x(1 : n+1);
% tangential forces
Tavector = x(n+2 : 2*(n+1));
% normal forces
Navector = x(2*(n+1)+1 : 3*(n+1));
% global orientations
psi = zeros(1,n);
for j = 1:n
    psi(j) = sum(phivector(1:j));
end

% iterate the segments
for i = 1:n
    % get the unknowns from vector x
    phi = phivector(i);
    Ta = Tavector(i);
    Na = Navector(i);
    phi_next = phivector(i+1);
    Ta_next = Tavector(i+1);
    Na_next = Navector(i+1);
    % relaxed state angles
    relax = phi_r(i);
    relax_next = phi_r(i+1);
    % toruqes
    Ma = -(phi-relax)*kt(i);
    Mb = (phi_next-relax_next)*kt(i+1);
    % force transformation
    Tb = -(Ta_next*cos(phi_next) - Na_next*sin(phi_next));
    Nb = -(Ta_next*sin(phi_next) + Na_next*cos(phi_next));

    % global Y coordinate
    Y= sum(L.*cos(psi(1:i-1)));
    % Y coordinate of contact point
    Yc = Y + L/2*cos(psi(i)) + h*cos(psi(i)-pi/2);
    % distance from contact surface
    delta = (D/2-d)-Yc;
    % give an approximation for the contact force
    Fy = epsilon*exp(-delta/epsilon);
    % the global horizontal component Fx = f(Fy)
    Fx = mu*Fy*tanh(10*v(0));
    % conversion to local coordinates
    Cx = Fx*cos(psi(i)- pi/2) - Fy*sin(psi(i) -pi/2);
    Cy = Fx*sin(psi(i) -pi/2) + Fy*cos(psi(i) -pi/2);
    
    % equations with contact point
    eqx = Ta + Tb + Cx;
    eqy = Na + Nb + Cy;
    eqt = Ma + Mb + L*Nb + h*Cx + L/2*Cy;

    
    % save equations
    eq = [eq, eqx, eqy, eqt];
end

% relative angle constraint
eqangle = sum(phivector)-pi;
% global geometry constraints
eqX = sum(L.*sin(psi))-D;
eqY = sum(L.*cos(psi));
% save constraints
eq = [eq, eqangle, eqX, eqY];
end
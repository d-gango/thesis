clear all

global n D d L kt h
n = 10;
contacts = []; % list of segments with assumed contact
D = 40; % diameter
d = 0; % contact depth
L = d*sin(pi/(2*n)); % segment length
kt = 1; % spring stiffness
h = 0.1; % lenght of contact segment

% initial values
x0 = zeros(1,(n+1)*3 + length(contacts));

x = fsolve(@equations, x0);


function eq = equations(x)
global n D d L kt h
eq = [];
for i = 1:n
    index = (i-1)*3; % x vector jump in ponit 
    % get the unknowns from vector x
    phi = x(index+1);
    Ta = x(index+2);
    Na = x(index+3);
    phi_next = x(index+4);
    Ta_next = x(index+5);
    Na_next = x(index+6);
    % toruqes
    Ma = -phi*kt;
    Mb = phi_next*kt;
    % force transformation
    Tb = -(Ta_next*cos(phi_next) - Na_next*sin(phi_next));
    Nb = -(Ta_next*sin(phi_next) + Na_next*cos(phi_next));
    
    % equations
    eqx = Ta + Tb;
    eqy = Na + Nb;
    eqt = Ma + Mb + L*Nb;
    
    % save equations
    eq = [eq, eqx, eqy, eqt];
end
% relative angles
phi = x(1:3:3*n+1);
% global orientations
psi = zeros(1,n);
for i = 1:n
    psi(i) = sum(phi(1:n));
end
% relative angle constraint
eqangle = sum(phi)-pi;
% global geometry constraints
eqX = sum(L.*sin(psi))-D;
eqY = sum(L.*cos(psi));

eq = [eq, eqangle, eqX, eqY];
end

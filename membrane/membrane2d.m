clear all

global n D d L kt h phi_r
n = 100;
contacts = []; % list of segments with assumed contact
D = 40; % diameter
d = 0; % contact depth
L = D*sin(pi/(2*n)); % segment length
kt = 1; % spring stiffness
h = 0.1; % lenght of contact segment
% angles at relaxed state
ang = pi/n;
phi_r = ones(1,n+1)*ang;
phi_r(1) = ang/2;
phi_r(end) = ang/2;

% initial values
x0 = zeros(1,(n+1)*3 + length(contacts));
% %-----------------------------------------
% % use relaxed state as initial condition
% for i = 1:n+1
%     x0(3*(i-1)+1) = phi_r(i);
% end
% %----------------------------------------
% options = optimoptions('fsolve','MaxFunctionEvaluations',1000000,...
%                         'MaxIterations', 5000);
x = fsolve(@equations, x0);

sol = equations(x);
drawsensor(x);

function eq = equations(x)
global n D d L kt h phi_r
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
    % relaxed state angles
    relax = phi_r(i);
    relax_next = phi_r(i+1);
    % toruqes
    Ma = -(phi-relax)*kt;
    Mb = (phi_next-relax_next)*kt;
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
    psi(i) = sum(phi(1:i));
end
% relative angle constraint
eqangle = sum(phi)-pi;
% global geometry constraints
eqX = sum(L.*sin(psi))-D;
eqY = sum(L.*cos(psi));

eq = [eq, eqangle, eqX, eqY];
end

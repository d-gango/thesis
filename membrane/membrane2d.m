clear all

global n D d L kt h phi_r contacts
n = 10;
contacts = [5,6]; % list of segments with assumed contact
D = 40; % diameter
d = 0.2; % contact depth
L = D*sin(pi/(2*n)); % segment length
kt = 1; % spring stiffness
h = 1; % lenght of contact segment
% angles at relaxed state
ang = pi/n;
phi_r = ones(1,n+1)*ang;
phi_r(1) = ang/2;
phi_r(end) = ang/2;

% initial values
x0 = zeros(1,(n+1)*3 + length(contacts));
% %-----------------------------------------
% use relaxed state as initial condition
for i = 1:n+1
    x0(3*(i-1)+1) = phi_r(i);
end
% %----------------------------------------
options = optimoptions('fsolve','MaxFunctionEvaluations',80000,...
                        'MaxIterations', 5000);
x = fsolve(@equations, x0, options);

sol = equations(x);
drawsensor(x);

function eq = equations(x)
global n D d L kt h phi_r contacts
eq = [];
% relative angles
phivector = x(1:3:3*n+1);
% global orientations
psi = zeros(1,n);
% iterate the segments
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
    
    for j = 1:n
        psi(j) = sum(phivector(1:j));
    end
    
    if ismember(i, contacts)
        [one, ci] = ismember(i, contacts);
        % the unknown global vertical contact force
        Fy = x(3*(n+1)+ci);
        % the global gorizontal component Fx = f(Fy)
        Fx = 0;
        % conversion to local coordinates
        Cx = Fx*cos(pi/2 - psi(i)) - Fy*sin(pi - psi(i));
        Cy = Fx*sin(pi/2 - psi(i)) + Fy*cos(pi - psi(i));
        
        % equations with contact point
        eqx = Ta + Tb + Cx;
        eqy = Na + Nb + Cy;
        eqt = Ma + Mb + L*Nb + h*Cx + L/2*Cy;
        
    else
        % equations with no cotact poit
        eqx = Ta + Tb;
        eqy = Na + Nb;
        eqt = Ma + Mb + L*Nb;
    end
    
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

for k = contacts
    % global Y coordinate
    Yk= sum(L.*cos(psi(1:k-1)));
    % Y coordinate of contact point
    Yck = Yk + L/2*cos(psi(k)) + h*cos(psi(k)-pi/2);
    % the constraint equation
    eqcont = D/2 - d - Yck;
    % save constraint
    eq = [eq, eqcont];
end
end

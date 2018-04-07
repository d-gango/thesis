clear all

n = 20;
contacts = [10 11];
% unknowns
phivector = sym('phi',[n+1 1]);
Tavector = sym('Ta', [n+1 1]);
Navector = sym('Na', [n+1 1]);
Fyvector = sym('Fy', [length(contacts) 1]);
% make them time dependent
syms t
for i = 1:n+1
    index = int2str(i);
    phiname = ['phi' index '(t)'];
    Taname = ['Ta' index '(t)'];
    Naname = ['Na' index '(t)'];
    phivector(i) = str2sym(phiname);
    Tavector(i) = str2sym(Taname);
    Navector(i) = str2sym(Naname);
end
for  i = 1:length(contacts)
    index = int2str(contacts(i));
    Fyname = ['Fy' index '(t)'];
    Fyvector(i) = str2sym(Fyname);
end
% constant variables
syms D d L kt h mu theta v b
phi_r = sym('phi_r', [1 n+1]);
% save the state variables
vars = [phivector; Tavector; Navector; Fyvector];

% global orientation
psi = sym('psi', [n 1]);
for i = 1:n
    psi(i) = sum(phivector(1:i));
end

eqns = [];
% iterate the segments
for i = 1:n
    % set the variables
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
    Ma = -(phi-relax)*kt - b*diff(phi);
    Mb = (phi_next-relax_next)*kt + b*diff(phi_next);
    % force transformation
    Tb = -(Ta_next*cos(phi_next) - Na_next*sin(phi_next));
    Nb = -(Ta_next*sin(phi_next) + Na_next*cos(phi_next));
    
    if ismember(i, contacts)
        [one, ci] = ismember(i, contacts);
        % the unknown global vertical contact force
        Fy = Fyvector(ci);
        % the global gorizontal component Fx = f(Fy)
        Fx = mu*Fy*tanh(100*v);
        % conversion to local coordinates
        Cx = Fx*cos(psi(i)- pi/2) - Fy*sin(psi(i) -pi/2);
        Cy = Fx*sin(psi(i) -pi/2) + Fy*cos(psi(i) -pi/2);
        
        % equations with contact point
        eqx = Ta + Tb + Cx == 0;
        eqy = Na + Nb + Cy == 0;
        eqt = Ma + Mb + L*Nb + h*Cx + L/2*Cy == theta*diff(phi,2);
        
    else
        % equations with no cotact poit
        eqx = Ta + Tb == 0;
        eqy = Na + Nb == 0;
        eqt = Ma + Mb + L*Nb == theta*diff(phi,2);
    end
    
    % save equations
    eqns = [eqns, eqx, eqy, eqt];
end
% relative angle constraint
eqangle = sum(phivector) == pi;
% global geometry constraints
eqX = sum(L.*sin(psi)) == D;
eqY = sum(L.*cos(psi)) == 0;
% save constraints
eqns = [eqns, eqangle, eqX, eqY];

for k = contacts
    % global Y coordinate
    Yk= sum(L.*cos(psi(1:k-1)));
    % Y coordinate of contact point
    Yck = Yk + L/2*cos(psi(k)) + h*cos(psi(k)-pi/2);
    % the constraint equation
    eqcont = D/2 - d == Yck;
    % save constraint
    eqns = [eqns, eqcont];
end

% add angular velocities to the variables
% add the corresponding equations
[eqns, vars, R] = reduceDifferentialOrder(eqns, vars);
%Reduce the differential index of the DAEs described by eqns and vars
[DAEs,DAEvars] = reduceDAEIndex(eqns,vars);
%Often, reduceDAEIndex introduces redundant equations and variables
%that can be eliminated. Eliminate redundant equations and variables
[DAEs,DAEvars] = reduceRedundancies(DAEs,DAEvars);
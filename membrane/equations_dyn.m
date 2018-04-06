function eq = equations_dyn(t,x)
global n D d L kt b h phi_r contacts mu theta v
eq = [];
% relative angles
phivector = x(1:n+1);
phidotvector = x(n+2:2*(n+1));
Tavector = x(2*(n+1)+1:3*(n+1));
Navector = x(3*(n+1)+1:4*(n+1));
% global orientations
psi = zeros(1,n);

% phi' = phidot
for i = 1:n+1
    % variables
    phidot = phidotvector(i);
    
    eqvel = phidot;
    eq = [eq; eqvel];
end

% phidot' = ...
for i = 1:n
    % variables
    phi = phivector(i);
    phidot = phidotvector(i);
    phi_next = phivector(i+1);
    phidot_next = phidotvector(i+1);
    Ta_next = Tavector(i+1);
    Na_next = Navector(i+1);
    % relaxed state angles
    relax = phi_r(i);
    relax_next = phi_r(i+1);
    % toruqes
    Ma = -(phi-relax)*kt - phidot*b;
    Mb = (phi_next-relax_next)*kt + phidot_next*b;
    % force transformation
    Nb = -(Ta_next*sin(phi_next) + Na_next*cos(phi_next));
    
    for j = 1:n
        psi(j) = sum(phivector(1:j));
    end
    
    if ismember(i, contacts)
        [one, ci] = ismember(i, contacts);
        % the unknown global vertical contact force
        Fy = x(4*(n+1)+ci);
        % the global horizontal component Fx = f(Fy)
        Fx = mu*Fy*tanh(100*v(t));
        % conversion to local coordinates
        Cx = Fx*cos(psi(i)- pi/2) - Fy*sin(psi(i) -pi/2);
        Cy = Fx*sin(psi(i) -pi/2) + Fy*cos(psi(i) -pi/2);
        
        % equations with contact point
        eqz = (Ma + Mb + L*Nb + h*Cx + L/2*Cy)/theta;
    else
        % equations with no cotact poit
        eqz = (Ma + Mb + L*Nb)/theta;
    end
    eq = [eq; eqz];
end
% static x
for i = 1:n
    % variables
    Ta = Tavector(i);
    phi_next = phivector(i+1);
    Ta_next = Tavector(i+1);
    Na_next = Navector(i+1);
    % force transformation
    Tb = -(Ta_next*cos(phi_next) - Na_next*sin(phi_next));
    
    for j = 1:n
        psi(j) = sum(phivector(1:j));
    end
    
    if ismember(i, contacts)
        [one, ci] = ismember(i, contacts);
        % the unknown global vertical contact force
        Fy = x(4*(n+1)+ci);
        % the global horizontal component Fx = f(Fy)
        Fx = mu*Fy*tanh(100*v(t));
        % conversion to local coordinates
        Cx = Fx*cos(psi(i)- pi/2) - Fy*sin(psi(i) -pi/2);
        
        % equations with contact point
        eqx = Ta + Tb + Cx;
    else
        % equations with no cotact poit
        eqx = Ta + Tb;
    end
    eq = [eq; eqx];
end
% static y
for i = 1:n
    % variables
    Na = Navector(i);
    phi_next = phivector(i+1);
    Ta_next = Tavector(i+1);
    Na_next = Navector(i+1);
    % force transformation
    Nb = -(Ta_next*sin(phi_next) + Na_next*cos(phi_next));
    
    for j = 1:n
        psi(j) = sum(phivector(1:j));
    end
    
    if ismember(i, contacts)
        [one, ci] = ismember(i, contacts);
        % the unknown global vertical contact force
        Fy = x(4*(n+1)+ci);
        % the global horizontal component Fx = f(Fy)
        Fx = mu*Fy*tanh(100*v(t));
        % conversion to local coordinates
        Cy = Fx*sin(psi(i) -pi/2) + Fy*cos(psi(i) -pi/2);
        
        % equations with contact point
        eqy = Na + Nb + Cy;
    else
        % equations with no cotact poit
        eqy = Na + Nb;
    end
    eq = [eq; eqy];
end

% relative angle constraint
eqangle = sum(phivector)-pi;
% global geometry constraints
eqX = sum(L.*sin(psi))-D;
eqY = sum(L.*cos(psi));
% save constraints
eq = [eq; eqangle; eqX; eqY];

for k = contacts
    % global Y coordinate
    Yk= sum(L.*cos(psi(1:k-1)));
    % Y coordinate of contact point
    Yck = Yk + L/2*cos(psi(k)) + h*cos(psi(k)-pi/2);
    % the constraint equation
    eqcont = D/2 - d - Yck;
    % save constraint
    eq = [eq; eqcont];
end
end
function eq = equations_slip_approx(x)
global n D d L kt h phi_r contacts mu eps
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
    
    % global Y coordinate
    Y= sum(L.*cos(psi(1:i-1)));
    % Y coordinate of contact point
    Yc = Y + L/2*cos(psi(i)) + h*cos(psi(i)-pi/2);
    % distance from contact surface
    delta = (D/2-d)-Yc;
    % give an approximation for the contact force
    Fy = eps*exp(-delta/eps);
    % the global horizontal component Fx = f(Fy)
    Fx = mu*Fy;
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
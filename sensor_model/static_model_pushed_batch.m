clear all
par = param();
filename = ['static_', num2str(par.n)];
load(filename);
% set grid for d
depth = 6;
F = [];
successful = [];

global contacts
par = param();
n = par.n;
phi_r = par.phi_r;
D = par.D;
if par.batch ~= 1
    error('Not in batch mode!')
end
if par.v(0) ~= 0
    error('Velocity should be 0!')
end


for d = depth
    % set global variable
    setGlobald(d);
    % initial values
    x0 = zeros(1,(n+1)*3);
    % %-----------------------------------------
    % use relaxed state as initial condition
    x0(1:n+1) = phi_r;
    % %----------------------------------------
    options = optimoptions('fsolve','MaxFunctionEvaluations',80000,...
        'MaxIterations', 5000);
    [x, fval, exitflag] = fsolve(@static_equations_approx, x0, options);
    if exitflag < 1
        error('Approximate solution falied')
    end
    
    % figuring out the contact points =========================================
    % approximate phi angles
    phi_a = x(1:n+1);
    % approximate orientations
    psi_a = getPsi(phi_a);
    % vertical positions of contact points
    Yc_a = getYc(psi_a);
    % guess the contact points
    ma = max(Yc_a);
    contact_bw = 0.25;
    mi = ma - contact_bw;
    % first guess
    contact_guess = find(Yc_a <= ma & Yc_a >= mi);
    % assume it's symmetric, take the first half
    contact_guess_half = contact_guess(1:length(contact_guess)/2);
    
    % try to find soultion with a subset of the guessed contact points
    cn = 0;
    solved = 0;
    while cn <= length(contact_guess_half)
        % get all combinations
        combinations = combnk(contact_guess_half, cn);
        % mirror
        other_half = sort((n+1)-combinations, 2);
        poss_cont = [combinations, other_half];
        % iterate the possible contacts
        rows = size(poss_cont,1);
        for i = 1:rows
            contacts = poss_cont(i,:);
            x0 = zeros(1,(n+1)*3 + length(contacts)); % new vector of unknowns
            % use relaxed state as initial condition
            x0(1:n+1) = phi_r;
            
            % solve the system
            [x, fval, exitflag] = fsolve(@static_equations, x0, options);
            disp(['d = ' num2str(d)]);
            disp(['Contacts: ', num2str(contacts)]);
            if exitflag > 0
                disp('Solution converged.');
                
                % check overlap
                phi = x(1:n+1);
                psi = getPsi(phi);
                Yc = getYc(psi);
                overlap = find(round(Yc,2) > (D/2-d));
                if isempty(overlap)
                    disp('No overlap.')
                    % check negative contact forces
                    Fy = getFy(x);
                    negFy = find(Fy < 0);
                    if isempty(negFy)
                        disp('No negative contact forces.')
                        disp('SOLUTION FOUND!')
                        solved = 1;
                        break;
                    else
                        disp('Negative contact forces');
                    end
                else
                    disp('Overlap.')
                end
            else
                disp('No solution');
            end
        end
        if solved
            break;
        else
            cn = cn+1;
        end
    end
    
    ind = find(round([solutions.d],2) == d);
    if isempty(ind)
        ind = length([solutions.d])+1;
    end
    solutions(ind).sol = x;
    solutions(ind).F = sum(Fy);
    solutions(ind).success = solved;
    solutions(ind).d = d;
    solutions(ind).contacts = contacts;
    
    
    save(filename,'solutions');
end

save(filename,'par', '-append');
%=========================================================================
function psi = getPsi(phi)
n = length(phi);
psi = zeros(1,n);
for j = 1:n
    psi(j) = sum(phi(1:j));
end
end
%=========================================================================
function Yc = getYc(psi)
par = param();
n = par.n;
L = par.L;
h = par.h;
Yc = zeros(1,n);
for i = 1:n
    % joint Y coordinate
    Y = sum(L.*cos(psi(1:i-1)));
    % contact point
    Yc(i) = Y + L/2*cos(psi(i)) + h*cos(psi(i)-pi/2);
end
end
%=========================================================================
function Fy = getFy(sol)
global contacts
par = param();
n = par.n;
Fy = [];
if ~isempty(contacts)
    Fy = sol(3*(n+1)+1:end);
end
end

clear all
filename = 'd_F_20.mat';
load(filename);
global n D d L kt h phi_r contacts mu
n = params.n;
mu = params.mu;
D = params.D; % diameter
depth = 6;%5.5 : 0.1 : 6;
L = D*sin(pi/(2*n)); % segment length
kt = params.kt; % spring stiffness
h = params.h; % lenght of contact segment
% angles at relaxed state
ang = pi/n;
phi_r = ones(1,n+1)*ang;
phi_r(1) = ang/2;
phi_r(end) = ang/2;

% load stored results
for d = depth
index = find(round([solutions.d],2) == round(d,2));
contacts = solutions(index).contacts;
% use push down solution as initial condition
x0_init = solutions(index).sol;

% figuring out the contact points =========================================
% approximate phi angles
phi_a = getPhi(x0_init);
% approximate orientations
psi_a = getPsi(phi_a);
% vertical positions of contact points
Yc_a = getYc(psi_a);
% guess the contact points
ma = max(Yc_a);
contact_bw = 1;
mi = ma - contact_bw;
% first guess
contact_guess = find(Yc_a <= ma & Yc_a >= mi);

% try to find soultion with a subset of the guessed contact points
cn = 0;
solved = 0;
while cn <= length(contact_guess)
    % get all combinations
    combinations = combnk(contact_guess, cn);
    % iterate the possible contacts
    rows = size(combinations,1);
    for i = 1:rows
        contacts = combinations(i,:);
        x0 = zeros(1,(n+1)*3 + length(contacts)); % new vector of unknowns
        x0(1:(n+1)*3) = x0_init(1:(n+1)*3); % use the initial state
        % solve the system
        [x, fval, exitflag] = fsolve(@equations_slip, x0);
        disp(['d: ', num2str(d)]);
        disp(['Contacts: ', num2str(contacts)]);
        if exitflag > 0
            disp('Solution converged.');
            
            % check overlap
            phi = getPhi(x);
    
            psi = getPsi(phi);
            Yc = getYc(psi);
            overlap = find(round(Yc,3) > (D/2-d));
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

if solved
    solutions(index).slip_solution = x;
    solutions(index).slip_contacts = contacts;
    save(filename, 'solutions', 'params');
end
end


%=========================================================================
function psi = getPsi(phi)
global n
psi = zeros(1,n);
for j = 1:n
    psi(j) = sum(phi(1:j));
end
end
%=========================================================================
function phi = getPhi(sol)
global n
phi = zeros(1,n+1);
for i = 1:n+1
    index = (i-1)*3;
    phi(i) = sol(index+1);
end
end
%=========================================================================
function Yc = getYc(psi)
global n L h
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
global contacts n
Fy = [];
if ~isempty(contacts)
    Fy = sol(3*(n+1)+1:end);
end
end



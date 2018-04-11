clear all

global n D d L kt h phi_r contacts mu b v theta
n = 20;
mu = 0.5;
D = 40; % diameter
d = 3; % contact depth
L = D*sin(pi/(2*n)); % segment length
kt = 100; % spring stiffness
b = 10; % damping
m = 100/n; %mass of one segment
theta = 1/3*m*L^2; % moment of inertia
h = 0.3; % lenght of contact segment
v = @(t) t;
% angles at relaxed state
ang = pi/n;
phi_r = ones(1,n+1)*ang;
phi_r(1) = ang/2;
phi_r(end) = ang/2;

% load stored results
load d_F_20.mat
index = find(round([solutions.d],2) == d);
contacts = solutions(index).contacts;
% use push down solution as initial condition
x0_init = solutions(index).sol;
drawsensor(x0_init);
% title('initial state');
[eqfun, DAEvars] = equations_dyn_sym();

% set initial condition guess
x0 = rearrange(x0_init);
x0est = zeros(1, length(DAEvars));
x0est(1:length(x0)) = x0;          % x0 guess
xp0est = zeros(1, length(DAEvars)); % x'0 guess
% Create an option set that specifies numerical tolerances
% for the numerical search.
opt = odeset('RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));
% Find consistent initial values for the variables and their derivatives
x0_fixed = zeros(1, length(DAEvars));
x0_fixed(1:length(x0)) = ones(1, length(x0));
[x0, xp0] = decic(eqfun, 0, x0est, [], xp0est, [], opt);

% solve DAE
tspan = [0 1];
[tSol,xSol] = ode15i(eqfun, tspan, x0, xp0, opt);

%=========================================================================
function x0 = rearrange(init)
% adds 0 angular velocities to the initial conditions
global n
phi = getPhi(init);
Ta = getTa(init);
Na = getNa(init);
rest = init(3*(n+1)+1:end);
x0 = [phi, Ta, Na, rest];
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
function Ta = getTa(sol)
global n
Ta = zeros(1,n+1);
for i = 1:n+1
    index = (i-1)*3;
    Ta(i) = sol(index+2);
end
end
%=========================================================================
function Na = getNa(sol)
global n
Na = zeros(1,n+1);
for i = 1:n+1
    index = (i-1)*3;
    Na(i) = sol(index+3);
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



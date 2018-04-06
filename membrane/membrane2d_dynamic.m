clear all

global n D d L kt h phi_r contacts mu b v theta
n = 20;
mu = 0.5;
D = 40; % diameter
d = 3; % contact depth
L = D*sin(pi/(2*n)); % segment length
kt = 100; % spring stiffness
b = 10; % damping
m = 0.1/n; %mass of one segment
theta = 1/3*m*L^2; % moment of inertia
h = 0.3; % lenght of contact segment
v = @(t) heaviside(t);
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
title('initial state');
% add the angular velocities to the initial conditions
x0 = setInit(x0_init);
% formulate the mass matrix
M = massMatrix();

options = odeset('Mass', M);
tspan = [0, 0.1];
[t,y] = ode23t(@equations_dyn,tspan,x0,options);


%=========================================================================
function M = massMatrix()
global n contacts
cn = length(contacts);
M = zeros(4*(n+1)+cn, 4*(n+1)+cn);
M(1:2*n+1, 1:2*n+1) = eye(2*n+1);
end
%=========================================================================
function x0 = setInit(init)
% adds 0 angular velocities to the initial conditions
global n
phi = getPhi(init);
Ta = getTa(init);
Na = getNa(init);
rest = init(3*(n+1)+1:end);
phidot = zeros(1,n+1);
x0 = [phi, phidot, Ta, Na, rest];
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



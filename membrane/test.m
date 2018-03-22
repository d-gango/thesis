clear all

global n D d L kt h phi_r contacts
n = 50;
D = 40; % diameter
d = 6; % contact depth
L = D*sin(pi/(2*n)); % segment length
kt = 100; % spring stiffness
h = 0.3; % lenght of contact segment
% angles at relaxed state
ang = pi/n;
phi_r = ones(1,n+1)*ang;
phi_r(1) = ang/2;
phi_r(end) = ang/2;

% initial values
%x0 = zeros(1,(n+1)*3 + length(contacts));
x0 = zeros(1,(n+1)*3);
% %-----------------------------------------
% use relaxed state as initial condition
for i = 1:n+1
    x0(3*(i-1)+1) = phi_r(i);
end
% %----------------------------------------
options = optimoptions('fsolve','MaxFunctionEvaluations',80000,...
                        'MaxIterations', 5000);
[x, fval, exitflag] = fsolve(@equations_approx, x0, options);

drawsensor(x);

% sort out the soultuons
phisol = zeros(1,n+1);
Tasol = zeros(1,n+1);
Nasol = zeros(1,n+1);
for i = 1:n+1
    index = (i-1)*3;
    phisol(i) = x(index+1);
    Tasol(i) = x(index+2);
    Nasol(i) = x(index+3);
end
Fysol = [];
if ~isempty(contacts)
    Fysol = x(3*(n+1)+1:end);
end



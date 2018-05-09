clear all
global contacts
contacts = [7 8 13 14]; % list of segments with assumed contact
par = param();
n = par.n;
phi_r = par.phi_r;
if par.batch ~= 0
    error('Should not be is batch mode!')
end
% initial values
x0 = zeros(1,(n+1)*3 + length(contacts));
% %-----------------------------------------
% use relaxed state as initial condition
x0(1:n+1) = phi_r;
% %----------------------------------------
options = optimoptions('fsolve','MaxFunctionEvaluations',80000,...
                        'MaxIterations', 5000, 'FunctionTolerance', 1e-8);
[x, fval, exitflag] = fsolve(@static_equations, x0, options);
if exitflag < 1
    disp('Solution failed')
end
animateSensor(0,x);

% sort out the soultuons
phisol = x(1 : n+1);
Tasol = x(n+2 : 2*(n+1));
Nasol = x(2*(n+1)+1 : 3*(n+1));
Fysol = [];
if ~isempty(contacts)
    Fysol = x(3*(n+1)+1:end);
end



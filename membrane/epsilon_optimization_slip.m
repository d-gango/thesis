clear all
filename = 'd_F_20.mat';
load(filename);
depth = cell2mat({solutions.d});
% initialisation
global n D d L kt h phi_r contacts mu
n = params.n;
D = params.D; % diameter
L = D*sin(pi/(2*n)); % segment length
kt = params.kt; % spring stiffness
h = params.h; % lenght of contact segment
mu = params.mu;
% angles at relaxed state
ang = pi/n;
phi_r = ones(1,n+1)*ang;
phi_r(1) = ang/2;
phi_r(end) = ang/2;
% use relaxed state as initial condition
x0 = zeros(1,(n+1)*3);
for i = 1:n+1
    x0(3*(i-1)+1) = phi_r(i);
end

% load solution from saved data
for d = 3
    index = find(round([solutions.d],2) == round(d,2));
    sol = solutions(index).slip_solution;
    contacts = solutions(index).slip_contacts;
    phi_sol = getPhi(sol);
    
    % iterate epsilon
    % initial interval
    start = 0;
    finish = 2;
    % iterate for dec decimal values
    for dec = 1:4
        res = 1/(10^dec);
        if round(start,5) == 0
            start = start+res;
        end
        epsvector = start : res : finish;
        errorvector = zeros(1,length(epsvector));
        for i = 1:length(epsvector)
            global eps
            eps = epsvector(i);
            
            options = optimoptions('fsolve','MaxFunctionEvaluations',80000,...
                'MaxIterations', 5000);
            [x, fval, exitflag] = fsolve(@equations_slip_approx, x0, options);
            approximations(i,:) = x;
            drawsensor(x);
            title(num2str(eps));
            if exitflag > 0
                disp(['d: ' num2str(d)]);
                disp(['epsilon: ' num2str(eps)]);
                disp('Solution found');
                
                phi_appr = getPhi(x);
                errorvector(i) = sqrt(sum((phi_appr - phi_sol).^2));
            else
                errorvector(i) = inf;
                disp(['d: ' num2str(d)]);
                disp(['epsilon: ' num2str(eps)]);
                disp('No solution.');
            end
        end
        % find minimum
        [m,ind] = min(errorvector);
        % save best values
        bestepsilon = epsvector(ind);
        bestapprox = approximations(ind,:);
        besterror = m;
        % generate new interval
        start = bestepsilon-res;
        finish = bestepsilon+res;
    end
    %     drawsensor(sol);
    %     title('solution');
    %     drawsensor(bestapprox);
    %     title(num2str(bestepsilon));
    % save values
    solutions(index).epsilon = bestepsilon;
    solutions(index).approximation = bestapprox;
    solutions(index).approximation_error = besterror;
    save(filename, 'solutions', 'params');
    
end
%==========================================================================
function phi = getPhi(sol)
global n
phi = zeros(1,n+1);
for i = 1:n+1
    index = (i-1)*3;
    phi(i) = sol(index+1);
end
end
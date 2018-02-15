clear all;
close all;

% parameters
global m b k mu g finish n
n = 5; % number of masses
m = ones(1, n) * 1;
b = ones(1, n) * 1;
b(end) = 0;
k = ones(1, n) * 1;
k(end) = 0;
mu = 0.5;
g = 1;

start = 0;
finish = 10;
dt = 0.001;
x0 = zeros(2*n, 1);
options = odeset('Events',@events);

tout = start;
xout = x0.';
teout = [];
xeout = [];
ieout = [];
frout = [];
slip = ones(1,n) * -1;

fr = zeros(1,n);
while tout(end) < finish
    % set friction force
    for i = 1:n
        if abs(xout(end,n+i)) < 1e-10  % zero velocity
            fh = fhat(xout(end,:),tout(end));
            if abs(fh(i)) > mu*m(i)*g
                slip(i) = 1;
                fr(i) = mu*m(i)*g*sign(fh(1));
            else
                slip(i) = 0;
                fr(i) = fh(i);
            end

        else
            slip(i) = 1;
            fr(i) = mu*m(i)*g*sign(xout(end,i+n));
        end
    end
    
    if all(slip)
   % solve until v=0
       [t,x,te,xe,ie] = ode45(@(t,x) eqOfMotion(x,t,fr),[start finish],...
                        x0, options);
    else
        slipping = find(slip);
        stuck = find(~slip);
        
        if isempty(slipping)
            t = [tout(end); tout(end)+dt];
            x = [xout(end,:); xout(end,:)];
            % no event
            te = [];
            xe = [];
            ie = [];
        else
            % solve the ode because of the slipping, but only for 1 step
            [t,x,te,xe,ie] = ode45(@(t,x) eqOfMotion(x,t,fr),[start start+dt],...
                        x0, options);
                    
            t = [tout(end); tout(end)+dt];
            x = [xout(end,:); x(end,:)];
            x(end,stuck) = xout(end, stuck);
            x(end,stuck+n) = xout(end, stuck+n);
%             if ~isempty(te)
%                 error('Another mass got stuck while solving the ode')
%             end
            te = [];
            xe = [];
            ie = [];
        end
    end

   % write output
   tout = [tout; t(2:end)];
   xout = [xout; x(2:end,:)];
   teout = [teout; te];
   xeout = [xeout; xe];
   ieout = [ieout; ie];
   if isempty(frout)
       frout = ones(size(t)) * fr;
   else
       frout = [frout; ones(size(t(2:end)))*fr];
   end
   
   
   % set the new ICs
   x0 = x(end,:)';
   
   start = t(end);
end

massnumbers = int2str((1:n)'); % for the labels

figure
plot(tout, xout(:,1:n));
xlim([0, finish]);
legend(massnumbers);
xlabel('t');
title('positions')

figure
plot(tout, xout(:,n+1:2*n));
xlim([0, finish]);
legend(massnumbers);
xlabel('t');
title('velocities')

figure
plot(xout(:,1:n), xout(:,n+1:2*n));
legend(massnumbers);
title('phase portraits')

figure
plot(tout, frout)
legend(massnumbers);
title('friction forces')
% --------------------------------------------------------------------------

function xdot = eqOfMotion(x,t,fr)
    global m n
    fh = fhat(x,t);
    xdot = zeros(2*n,1);
    for i = 1:n
        xdot(i) = x(n+i);
        xdot(n+i) = 1/m(i)*(fh(i)-fr(i));
    end
end

% --------------------------------------------------------------------------

function fh = fhat(x, t)
    global b k n
    fh = zeros(n,1);
    fh(1) = forceIn(t) - k(1)*(x(1)-x(2)) - b(1)*(x(n+1)-x(n+2));
    if n>2
        for j = 2:n-1
            fh(j) = k(j-1)*(x(j-1)-x(j)) + b(j-1)*(x(n+j+1)-x(n+j))...
                   -k(j)*(x(j)-x(j+1)) - b(j)*(x(n+j)-x(n+j+1));
        end
    end
    fh(n) = k(n-1)*(x(n-1)-x(n)) + b(n-1)*(x(2*n-1)-x(2*n))...
           -k(n)*x(n) - b(n)*x(2*n);
end

% --------------------------------------------------------------------------

function f = forceIn(t)
    f = t;
end

% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,x)
    global n
    % when one velocity is zero stop the solver
    value = zeros(n,1);
    for i = 1:n
        value(n) = x(n+i);  % velocity = 0
    end
    isterminal = ones(n,1); % stop the integration
    direction = zeros(n,1);   % both directions
end
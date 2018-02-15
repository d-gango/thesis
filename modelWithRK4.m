clear all;
close all;

% parameters
global m b k mu g n
n = 10; % number of masses
m = ones(1, n) * 1;
b = ones(1, n) * 10;
b(end) = 0;
k = ones(1, n) * 1;
k(end) = 0;
mu = 0.5;
g = 1;

t = 0;
dt = 0.001;
finish = 15;

x0 = zeros(1, 2*n);
% calculate external and friction forces
fh = fhat(x0,t);
[fr, slip] = friction(x0, fh);

tout = t;
xout = x0;
fhout = fh;
frout = fr;
slipout = slip;

while t < finish
    % integrate 1 timestep
    x = RK4(@eqOfMotion, x0, t, dt);
    
    % set the new ICs
    x0 = x;
    t = t + dt;
    
    % calculate external and friction forces
    fh = fhat(x0,t);
    [fr, slip] = friction(x0, fh);

   % write output
   tout = [tout; t];
   xout = [xout; x];
   fhout = [fhout; fh];
   frout = [frout; fr];
   slipout = [slipout; slip];
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

function xdot = eqOfMotion(x,t)
    global m n
    fh = fhat(x,t);
    [fr, slip] = friction(x, fh);
    slipping = find(slip);
    xdot = zeros(1,2*n);
    if isempty(slipping)
        return;
    else
        for i = slipping
            xdot(i) = x(n+i);
            xdot(n+i) = 1/m(i)*(fh(i)-fr(i));
        end
    end
end

% --------------------------------------------------------------------------

function fh = fhat(x, t)
    global b k n
    fh = zeros(1,n);
    fh(1) = forceIn(t) - k(1)*(x(1)-x(2)) - b(1)*(x(n+1)-x(n+2));
    if n>2
        for j = 2:n-1
            fh(j) = k(j-1)*(x(j-1)-x(j)) + b(j-1)*(x(n+j-1)-x(n+j))...
                   -k(j)*(x(j)-x(j+1)) - b(j)*(x(n+j)-x(n+j+1));
        end
    end
    fh(n) = k(n-1)*(x(n-1)-x(n)) + b(n-1)*(x(2*n-1)-x(2*n))...
           -k(n)*x(n) - b(n)*x(2*n);
end

% --------------------------------------------------------------------------

function [fr, slip] = friction(x, fh)
    global n mu m g
    slip = ones(1,n) * -1;
    fr = zeros(1,n);
    for i = 1:n
        if abs(x(n+i)) < 1e-3  % zero velocity
            if abs(fh(i)) > mu*m(i)*g
                slip(i) = 1;
                fr(i) = mu*m(i)*g*sign(fh(i));
            else
                slip(i) = 0;
                fr(i) = fh(i);
            end

        else
            slip(i) = 1;
            fr(i) = mu*m(i)*g*sign(x(i+n));
        end
    end
end

% --------------------------------------------------------------------------

function f = forceIn(t)
    f = t;
end

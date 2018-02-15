clear all;
% parameters
global m b k mu g finish
m = 200;
b = 5;
k = 5000;
mu = 1;
g = 1;

start = 0;
finish = 2.5;
dt = 0.001;
x0 = [0; 0];
options = odeset('Events',@events);

tout = start;
xout = x0.';
teout = [];
xeout = [];
ieout = [];
frout = [];
slip = -1;

while tout(end) < finish
    % set friction force
    if abs(xout(end,2)) < 1e-6  % zero velocity
        if abs(fhat(xout(end,:),tout(end))) > mu*m*g
            slip = 1;
            fr = mu*m*g*sign(fhat(xout(end,:),tout(end)));
        else
            slip = 0;
            fr = fhat(xout(end,:),tout(end));
        end
        
    else
        slip = 1;
        fr = mu*m*g*sign(xout(end,2));
    end
    
    if slip
   % solve until v=0
       [t,x,te,xe,ie] = ode45(@(t,x) eqOfMotion(x,t,fr),[start finish],...
                        x0, options);
%        % A good guess of a valid first timestep is the length of the last valid
%        % timestep, so use it for faster computation.  'refine' is 4 by default.
%        options = odeset(options,'InitialStep',t(end)-t(end-refine),...
%       'MaxStep',t(end)-t(1));
    else
        t = [tout(end); tout(end)+dt];
        x = [xout(end,:); xout(end,:)];
        % no event
        te = [];
        xe = [];
        ie = [];
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
   x0(1) = x(end,1);
   x0(2) = x(end,2);
   
   start = t(end);
end

figure
plot(tout, xout(:,1));
xlim([0, finish]);
%ylim([-2, 2]);
xlabel('t');
ylabel('x');
title('position')

figure
plot(tout, xout(:,2));
xlim([0, finish]);
%ylim([-2, 2]);
xlabel('t');
ylabel('xdot');
title('velocity')

figure
plot(xout(:,1), xout(:,2));
% xlim([-0.001 0.0014])
% ylim([-0.03 0.03])
title('phase portrait')

figure
plot(xout(:,1), frout,'.')
% xlim([-0.0008 0.0013])
% ylim([-250 250])
title('friction force')
% --------------------------------------------------------------------------

function xdot = eqOfMotion(x,t,fr)
global m
xdot = zeros(2,1);
xdot(1) = x(2);
xdot(2) = 1/m*(fhat(x,t) - fr);
end

% --------------------------------------------------------------------------

function fh = fhat(x, t)
    global b k
    fh = forceIn(t) - b*x(2) - k*x(1);
end

% --------------------------------------------------------------------------

function f = forceIn(t)
    f = 300*sin(8*pi*t);
end

% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,x)
% when the velocity is zero or when the simulation time is over
% stop the solver
global finish
value = [x(2); finish-t];     % velocity = 0; t = finish
isterminal = [1; 1];   % stop the integration
direction = [0; 0];   % both directions
end
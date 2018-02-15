clear all;
% parameters
global m b k mu g finish
m = [1, 1];
b = [1, 1];
k = [1, 1];
mu = 0.5;
g = 1;

start = 0;
finish = 20;
dt = 0.001;
x0 = [0; 0; 0; 0];
options = odeset('Events',@events);

tout = start;
xout = x0.';
teout = [];
xeout = [];
ieout = [];
frout = [];
slip = [-1, -1];

fr = zeros(1,2);
while tout(end) < finish
    % set friction force
    for i = 1:2
        if abs(xout(end,i+2)) < 1e-10  % zero velocity
            fh = fhat(xout(end,:),tout(end));
            if abs(fh(i)) > mu*m(i)*g
                slip(i) = 1;
                fr(i) = mu*m(i)*g*sign(fh(i));
            else
                slip(i) = 0;
                fr(i) = fh(i);
            end

        else
            slip(i) = 1;
            fr(i) = mu*m(i)*g*sign(xout(end,i+2));
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
            x(end,stuck+2) = xout(end, stuck+2);
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

plot(tout, xout(:,1),'b', tout, xout(:,2),'r',...
     tout, xout(:,3),':b', tout, xout(:,4),':r');
xlim([0, finish]);
%ylim([-2, 2]);
% hold on;
% % plot events 
% if ~isempty(teout)
%     plot(teout,xeout(:,1),'ko')
% end
legend('x1', 'x2', 'x3', 'x4');
xlabel('t');
% hold off

figure
plot(xout(:,1), xout(:,3), 'b', xout(:,2), xout(:,4), 'r');
legend('1', '2');

figure
plot(tout, frout(:,1), 'b', tout, frout(:,2), 'r')
legend('1', '2');
title('friction force')
% --------------------------------------------------------------------------

function xdot = eqOfMotion(x,t,fr)
global m
fh = fhat(x,t);
xdot = zeros(4,1);
xdot(1) = x(3);
xdot(2) = x(4);
xdot(3) = 1/m(1)*(fh(1) - fr(1));
xdot(4) = 1/m(2)*(fh(2) - fr(2));
end

% --------------------------------------------------------------------------

function fh = fhat(x, t)
    global b k
    fh = [forceIn(t) - k(1)*(x(1)-x(2)) - b(1)*(x(3)-x(4));...
         k(1)*(x(1)-x(2)) + b(1)*(x(3)-x(4)) - k(2)*x(2) - b(2)*x(4)];
end

% --------------------------------------------------------------------------

function f = forceIn(t)
    f = 2*sin(t);
end

% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,x)
% when the velocity is zero or when the simulation time is over
% stop the solver
value = [x(3); x(4)];     % velocity = 0; t = finish
isterminal = [1; 1];   % stop the integration
direction = [0; 0];   % both directions
end
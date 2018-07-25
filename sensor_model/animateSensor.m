function animateSensor(tSol, xSol)
par = param();
D = par.D;
L = par.L;
n = par.n;
if par.batch
    d = getGlobald();
else
    d = par.d(0);
end
h = par.h;
global contacts
% relative angles
phi = xSol(end,1:n+1);
% global orientations
psi = zeros(1,n);
for i = 1:n
    psi(i) = sum(phi(1:i));
end
% global X coordinates
X = zeros(1,n+1);
for i = 1:n
    X(i+1)= sum(L.*sin(psi(1:i)));
end
% global Y coordinates
Y = zeros(1,n+1);
for i = 1:n
    Y(i+1)= sum(L.*cos(psi(1:i)));
end
% coordinates for drawing the contact thingies
s.x = [0 0];
s.y = [0 0];
contactlines = repmat(s, [1 n]);
for i = 1:n
    % coordinates of the segment endpoints
    XA = X(i);
    XB = X(i+1);
    YA = -Y(i);
    YB = -Y(i+1);
    % starting point of contact line
    x1 = (XA+XB)/2;
    y1 = (YA+YB)/2;
    % endpoint of contact line
    x2 = x1 + h*cos(psi(i)-pi);
    y2 = y1 + h*sin(psi(i)-pi);
    % save coordinates
    contactlines(i).x = [x1 x2];
    contactlines(i).y = [y1 y2];
end

figure
line([0 D], [0 0], 'Color', 'k') % sensor base
hold on
line([-10 D+10], [-D/2-h+d -D/2-h+d], 'Color', 'k') % contact surface
sensorplot = plot(X,-Y, 'LineWidth', 2); % sensor
for i = 1:n
    if ismember(i,contacts)
        contactcolor = 'red';
    else
        contactcolor = 'blue';
    end
    contactplot(i) = line(contactlines(i).x, contactlines(i).y,...
                'Color', contactcolor, 'LineWidth', 2);
end
timetext = text(D, -D/2, num2str(tSol(1)));
xlim([-10, D+10]);
ylim([-D/2-5, 5]);
axis equal
%grid on

for j = 1:length(tSol)
    % relative angles
    phi = xSol(j,1:n+1);
    % global orientations
    for i = 1:n
        psi(i) = sum(phi(1:i));
    end
    % global X coordinates
    for i = 1:n
        X(i+1)= sum(L.*sin(psi(1:i)));
    end
    % global Y coordinates
    for i = 1:n
        Y(i+1)= sum(L.*cos(psi(1:i)));
    end
    for i = 1:n
        % coordinates of the segment endpoints
        XA = X(i);
        XB = X(i+1);
        YA = -Y(i);
        YB = -Y(i+1);
        % starting point of contact line
        x1 = (XA+XB)/2;
        y1 = (YA+YB)/2;
        % endpoint of contact line
        x2 = x1 + h*cos(psi(i)-pi);
        y2 = y1 + h*sin(psi(i)-pi);
        % save coordinates
        contactlines(i).x = [x1 x2];
        contactlines(i).y = [y1 y2];
    end
    set(sensorplot, 'XData', X, 'YData', -Y);
    for i = 1:n
        set(contactplot(i), 'XData', contactlines(i).x,...
                            'YData', contactlines(i).y);
    end
    set(timetext, 'String', ''); %set(timetext, 'String', num2str(tSol(j)));
    drawnow
    
end
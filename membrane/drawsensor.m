function drawsensor(x)
global D L n d contacts h
% relative angles
phi = x(1:3:3*n+1);
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
if d>0
line([-10 D+10], [-D/2+d -D/2+d], 'Color', 'k') % contact surface
end
plot(X,-Y) % sensor
for i = 1:n
    if ismember(i,contacts)
        contactcolor = 'red';
    else
        contactcolor = 'blue';
    end
    line(contactlines(i).x, contactlines(i).y, 'Color', contactcolor);
end
xlim([-10, D+10]);
ylim([-D/2-5, 5]);
axis equal
grid on


end
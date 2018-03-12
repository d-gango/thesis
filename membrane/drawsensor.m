function drawsensor(x)
global D L n d
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

figure
line([0 D], [0 0], 'Color', 'k') % sensor base
hold on
if d>0
line([-10 D+10], [-D/2+d -D/2+d], 'Color', 'k') % contact surface
end
plot(X,-Y)
xlim([-10, D+10]);
ylim([-D/2-5, 5]);
axis equal
grid on


end
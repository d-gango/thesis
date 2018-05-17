function drawSurface(joint_coordinates, pin_coordinates)
par = param();
x = joint_coordinates(1,:);
y = joint_coordinates(2,:);
z = joint_coordinates(3,:);

xlin = linspace(min(x), max(x), 50);
zlin = linspace(min(z), max(z), 50);
[X,Z] = meshgrid(xlin, zlin);

Y = griddata(x,-z,y,X,-Z,'natural');
mesh(X,-Z,Y)
%axis([-par.D/2 par.D/2 -par.D/2  par.D/2 -par.D/2 0]);
axis equal
xlabel('x'); ylabel('-z'); zlabel('y')
hold on
plot3(x, -z, y, '.', 'MarkerSize', 5);
plot3(pin_coordinates(1,:), -pin_coordinates(3,:), pin_coordinates(2,:),...
     '.r', 'MarkerSize', 10);
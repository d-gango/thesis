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
% figure
hold on
scatter3(x, -z, y, '.');
% hold on
scatter3(pin_coordinates(1,:), -pin_coordinates(3,:), pin_coordinates(2,:),...
     'MarkerFaceColor', 'flat');
 
axis equal
zlim([-20 0]);
xlabel('x [mm]'); ylabel('-z [mm]'); zlabel('y [mm]')
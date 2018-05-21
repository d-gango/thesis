clear all

pins = getPinData();
load('tactip-benchmarks\trainDepthVideo05012246/expt.mat')
load('tactip-benchmarks\trainDepthVideo05012246/expt_mapping.mat')
load('sol.mat')
cam_reference = Expt.pinPositions';
cam_reference = cam_reference(:,expt_mapping);

figure
plot(cam_reference(1,:),cam_reference(2,:), '.')
axis equal
xlabel('x'); ylabel('y')
title('original camera image')

par = param();
[deformed_joints, deformed_pins] = deformedShape3D(phisol);
% figure
% plot3(deformed_pins(1,:),-deformed_pins(3,:),deformed_pins(2,:), '.')
% axis equal
% xlabel('x'); ylabel('-z'); zlabel('y')
% title('deformed pins')

% % interpolation
% cam_deformed_x = griddata(pins.original_coordinates(1,:),...
%                  pins.original_coordinates(2,:),...
%                  pins.original_coordinates(3,:),...
%                  cam_reference(1,:),...
%                  deformed_pins(1,:),deformed_pins(2,:),deformed_pins(3,:),...
%                  'linear');
% cam_deformed_y = griddata(pins.original_coordinates(1,:),...
%                  pins.original_coordinates(2,:),...
%                  pins.original_coordinates(3,:),...
%                  cam_reference(2,:),...
%                  deformed_pins(1,:),deformed_pins(2,:),deformed_pins(3,:),...
%                  'linear');

% interpolation
fx = scatteredInterpolant(pins.original_coordinates(1,:)',...
                 pins.original_coordinates(2,:)',...
                 pins.original_coordinates(3,:)',...
                 cam_reference(1,:)', 'natural' ,'linear');
fy = scatteredInterpolant(pins.original_coordinates(1,:)',...
                 pins.original_coordinates(2,:)',...
                 pins.original_coordinates(3,:)',...
                 cam_reference(2,:)', 'natural' ,'linear');
cam_deformed_x = fx(deformed_pins(1,:), deformed_pins(2,:), deformed_pins(3,:));
cam_deformed_y = fy(deformed_pins(1,:), deformed_pins(2,:), deformed_pins(3,:));

figure
plot(cam_deformed_x, cam_deformed_y, '.')
axis equal
xlabel('x'); ylabel('y')
title('interpolated deformed image')
% hold on
% for i = 1:length(cam_deformed_x)
%     plot(cam_deformed_x(i), cam_deformed_y(i), '.r', 'MarkerSize', 10)
% end

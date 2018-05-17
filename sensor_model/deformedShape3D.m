function [deformed_joints, deformed_pins] = deformedShape3D(phisol)
pins = getPinData();
par = param();
% calculate joint coordinates
psi = getPsi(phisol);
joints = zeros(3,par.n+1);
joints(1,1) = -par.D/2;
for i = 2:par.n+1
    joints(1,i) = joints(1,i-1) + par.L*sin(psi(i-1));
    joints(2,i) = joints(2,i-1) - par.L*cos(psi(i-1));
end
% generate 3D joint coordinates
deformed_joints = [];
res = 50;
d_alpha = linspace(0, pi, res);
% rotate cross section with d_alpha resolution
for i = 1:res-1
    T = rotationVectorToMatrix(d_alpha(i) * [0 1 0]);
    deformed_joints = [deformed_joints, T*joints];
end
%drawSurface(deformed_joints)
    
deformed_pins = zeros(3,length(pins.alpha));
for i = 1:length(pins.alpha)
    ind = pins.segment_index(i); % corresponding segment index
%     figure
%     plot(pins.relaxed_joints(1,:), pins.relaxed_joints(2,:), '-o')
%     axis equal; hold on
%     scatter(pins.plane_coordinates(1), pins.plane_coordinates(2));
    T = rotationVectorToMatrix(-(phisol(ind)-par.phi_r(ind)) * [0  0 1]);
%     plot(joints(1,:), joints(2,:), '-o')
%     axis equal; hold on
    deformed_planar_coordinates = ...
        joints(:,ind) + T*[pins.position_vector(:,i);0];
%     scatter(deformed_planar_coordinates(1), deformed_planar_coordinates(2))
    R = rotationVectorToMatrix(pins.alpha(i) * [0 1 0]);
    deformed_pins(:,i) = R * deformed_planar_coordinates;
%     rotated_joints = R * joints;
%     figure
%     plot3(rotated_joints(1,:), -rotated_joints(3,:), rotated_joints(2,:),'-o')
%     axis equal; hold on;
%     scatter3(deformed_pins(1,i), -deformed_pins(3,i), deformed_pins(2,i))
%     xlabel('x'); ylabel('-z'); zlabel('y')
end
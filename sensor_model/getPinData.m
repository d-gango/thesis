% loads 3D pin coordinates from file pin_data.mat
% calculates the pins relative angle to the XY plane
% finds the closest segment of the 2D model to each pin
function pins = getPinData()
load pin_data.mat
points = round(points, 2);
pins.original_coordinates = points;
% original configuration
% figure
% scatter3(points(1,:), -points(3,:), points(2,:), 'MarkerFaceColor', 'flat');
% %title('Initial pin positions')
% axis equal
% zlim([-20 0])
% xlabel('x [mm]'); ylabel('-z [mm]'); zlabel('y [mm]')

par = param();
psi = getPsi(par.phi_r);
% coordinates of joints in the cross-section plane
joint_x = zeros(1, par.n+1);
joint_x(1) = -par.D/2; % shift the 0 of x axis to the left
joint_y = zeros(1, par.n+1);
for i = 2:par.n+1
    joint_x(i) = joint_x(i-1) + par.L*sin(psi(i-1));
    joint_y(i) = joint_y(i-1) - par.L*cos(psi(i-1));
end
pins.relaxed_joints = [joint_x; joint_y];
% find rotation angle of corresponding cross-section plane
pins.alpha = atan2(points(3,:), points(1,:));
for i = 1:length(pins.alpha)
    T = rotationVectorToMatrix(-pins.alpha(i)*[0 1 0]);
    p = T*points(:,i); % coordinates in cross section plane
    pins.plane_coordinates(:,i) = p(1:2);
    % calculate distance from joints
    distance = zeros(1, par.n+1);
    for j =  1:par.n+1
        distance(j) = norm([joint_x(j);joint_y(j)]-p(1:2)); 
    end
    % find the closest 2 joints
    increasing = sort(distance);
    joint1 = find(distance == increasing(1));
    joint2 = find(distance == increasing(2));
    if abs(joint1-joint2) ~= 1
        error('Not consecutive joints!')
    end
    pins.segment_index(i) = min([joint1, joint2]); % save corresponding segment index
    pins.position_vector(:,i) = p(1:2) - ...
          [joint_x(pins.segment_index(i)); joint_y(pins.segment_index(i))];
end
end
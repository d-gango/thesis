clear all
% load solution file
load('static_20_thick.mat');
% load measurement data
measurement_data = getMeasurementData();
% get interpolation functions
[fx, fy] = getInterpolatingFun();

% generate 2D data from simulation
for i = 1:length(measurement_data)
    % get 3D pin coordinates
    [deformed_joints, deformed_pins] = deformedShape3D(solutions(i).sol(1:21));
    simulation_data{i} = ...
        [fx(deformed_pins(1,:), deformed_pins(2,:), deformed_pins(3,:));...
        fy(deformed_pins(1,:), deformed_pins(2,:), deformed_pins(3,:))];
    % caluclate displacements
    measurement_displacement{i} = measurement_data{i} - measurement_data{1};
    simulation_displacement{i} = simulation_data{i} - simulation_data{1};
end

% plot every step
for i = 1:10:60
    figure
    plot(measurement_data{i}(1,:), measurement_data{i}(2,:), '.',...
         simulation_data{i}(1,:), simulation_data{i}(2,:), '.r',...
         measurement_data{1}(1,:), measurement_data{1}(2,:), '.k');
    axis equal
    xlabel('x'); ylabel('y');
end

% comparison
for i = 1:length(measurement_data)
    difference{i} = simulation_data{i} - measurement_data{i};
    abs_difference(:,i) = [sum(abs(difference{i}(1,:)));...
        sum(abs(difference{i}(2,:)))];
    displacement_rel_error{i} = ...
      (simulation_displacement{i} - measurement_displacement{i}) ./...
      measurement_displacement{i} * 100;
end

figure
plot(1:length(measurement_data), abs_difference(1,:),...
    1:length(measurement_data), abs_difference(2,:))
legend('x', 'y')
ylabel('error')

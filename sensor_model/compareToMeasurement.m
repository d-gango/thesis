function [quad_error, displacement_rel_error]= compareToMeasurement(sol,d)
% load measurement data
measurement_data = getMeasurementData();
par = param();
% get interpolation functions
[fx, fy] = getInterpolatingFun();

% get corresponding index in the measurements
index = d/0.1 +1;


% generate 2D data from simulation
for i = 1:length(d)
    % get 3D pin coordinates
    [deformed_joints, deformed_pins] = deformedShape3D(sol(i,1:par.n+1));
    simulation_data{i} = ...
        [fx(deformed_pins(1,:), deformed_pins(2,:), deformed_pins(3,:));...
        fy(deformed_pins(1,:), deformed_pins(2,:), deformed_pins(3,:))];
    % caluclate displacements and error
    measurement_displacement{i} = measurement_data{index(i)} - measurement_data{1};
    simulation_displacement{i} = simulation_data{i} - measurement_data{1};
    simulation_error{i} = simulation_data{i} - measurement_data{index(i)};
end

% plot every step
for i = 1:length(d)
    figure
    plot(measurement_data{index(i)}(1,:), measurement_data{index(i)}(2,:), '.b',...
        simulation_data{i}(1,:), simulation_data{i}(2,:), '.r',...
        measurement_data{1}(1,:), measurement_data{1}(2,:), '.k', 'MarkerSize', 10);
    axis equal
    xlabel('x'); ylabel('y');
    title(['d = ', num2str(d(i)), ' mm']);
    legend('measurement', 'simulation', 'original')
end

% comparison
for i = 1:length(d)
    quad_error{i} = sqrt(sum((simulation_data{i} - measurement_data{index(i)}).^2, 2))
    for j = 1:length(simulation_displacement{i})
        displacement_rel_error{i}(j) = ...
            (norm(simulation_displacement{i}(:,j)) - norm(measurement_displacement{i}(:,j))) ./...
            norm(measurement_displacement{i}(:,j)) * 100;
    end
    disp('quadratic error:');
    disp(quad_error{i});
end

% find pins with small relative errors
for i = 1:length(d)
    treshold = 25;
    closest = find(displacement_rel_error{i} < treshold & displacement_rel_error{i} > -treshold);
    figure
    plot(measurement_data{1}(1,:), measurement_data{1}(2,:), '.k',...
        measurement_data{1}(1,closest), measurement_data{1}(2,closest), '.r','MarkerSize', 12);
    axis equal
    xlabel('x'); ylabel('y');
    title(['d = ', num2str(d(i)), ' mm']);
    legend('all', 'good approximations')
end

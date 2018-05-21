pins = getPinData();
load('tactip-benchmarks\trainDepthVideo05012246/expt.mat')
load('tactip-benchmarks\trainDepthVideo05012246/expt_mapping.mat')
figure
plot(Expt.pinPositions(:,1),Expt.pinPositions(:,2), '.')
axis equal
xlabel('x'); ylabel('y')
title('original camera image')

hold on
for i = expt_mapping
    plot(Expt.pinPositions(i,1), Expt.pinPositions(i,2), '.r', 'MarkerSize', 15);
end

% figure
% scatter3(pins.original_coordinates(1,:), -pins.original_coordinates(3,:),...
%     pins.original_coordinates(2,:));
% axis equal
% xlabel('x'); ylabel('-z'); zlabel('y')

f = 150;
camera_position = 50;
camera_coordinates = pins.original_coordinates;
camera_coordinates(2,:) =camera_coordinates(2,:)-camera_position;

% figure
% scatter3(camera_coordinates(1,:), -camera_coordinates(3,:),...
%     camera_coordinates(2,:));
% axis equal
% xlabel('x'); ylabel('-z'); zlabel('y')

projected_coordinates =zeros(3,length(pins.alpha));
for i = 1:length(pins.alpha)
    projected_coordinates(:,i) = ...
        f/camera_coordinates(2,i) * camera_coordinates(:,i);
end
figure
plot(camera_coordinates(1,:),-camera_coordinates(3,:), '.')
axis equal
xlabel('x'); ylabel('z');
title('perspective projection')
hold on
for i = 1:length(camera_coordinates(1,:))
    plot(camera_coordinates(1,i), -camera_coordinates(3,i), '.r', 'MarkerSize', 15);
end
F = cell2mat({solutions.F});
depth = cell2mat({solutions.d});
success = cell2mat({solutions.success});

figure
plot(depth,F);
axis equal
xlabel('d [mm]');
ylabel('F [N]');

fitted_parab = fit(depth',F','poly2');
fitted_lin = fit(depth',F','poly1');

figure
plot(fitted_parab, depth, F);
% hold on
% plot(fitted_lin, 'g');
xlabel('d [mm]');
ylabel('F [N]');
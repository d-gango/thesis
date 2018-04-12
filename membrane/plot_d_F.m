F = cell2mat({solutions.F});
depth = cell2mat({solutions.d})+0.3;
success = cell2mat({solutions.success});
epsilon = cell2mat({solutions.epsilon});
appr_error = cell2mat({solutions.approximation_error});

% figure
% plot(depth,F);
% axis equal
% xlabel('d [mm]');
% ylabel('F [N]');

fitted_parab = fit(depth',F','poly2');
fitted_lin = fit(depth',F','poly1');

figure
plot(depth, F,'.');
% hold on
% plot(fitted_lin, 'g');
xlabel('d [mm]');
ylabel('F_p [mN]');

figure
plot(depth, epsilon)
hold on
plot(depth, appr_error);
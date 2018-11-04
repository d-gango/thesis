close all
clear
load spinodal_data_F.mat
load coulomb_data_F.mat

figure
plot(spinodal.t,spinodal.x(:,1), '-',...
     coulomb.t, coulomb.x(:,1), '-', 'LineWidth', 2, 'MarkerSize', 15)
title('position')
legend('gen. spinodal law', ['Coulomb' char(39) 's law'])
ylabel('x [m]'); xlabel('t [s]');

figure
plot(spinodal.t,spinodal.x(:,2), '-',...
     coulomb.t,coulomb.x(:,2), '-','LineWidth', 2, 'MarkerSize', 15)
title('velocity')
legend('gen. spinodal law', ['Coulomb' char(39) 's law'])
ylabel('v [m/s]'); xlabel('t [s]');

figure
plot(spinodal.t,spinodal.x(:,3), '-', 'LineWidth', 2, 'MarkerSize', 15)
title('"surface roughness"')
ylabel('$\phi$ [1]', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('t [s]', 'Interpreter', 'latex', 'FontSize', 14);

figure
mu = coulomb.mu*ones(1,length(coulomb.t));
plot(spinodal.t,spinodal.mu, '-',...
     coulomb.t,mu, '-','LineWidth', 2, 'MarkerSize', 15)
title('friction coeffitient')
legend('gen. spinodal law', ['Coulomb' char(39) 's law'])
ylabel('\mu [1]'); xlabel('t [s]');

% figure
% plot(spinodal.t,spinodal.v_rel, '-',...
%      coulomb.t,coulomb.v_rel, '-','LineWidth', 2, 'MarkerSize', 15)
% title('relative velocity')
% legend('gen. spinodal law', ['Coulomb' char(39) 's law'])
% ylabel('v_{rel} [m/s]'); xlabel('t [s]');

figure
plot(spinodal.t,spinodal.Ff, '-',...
     coulomb.t,coulomb.Ff, '-', 'LineWidth', 2, 'MarkerSize', 15)
title('friction force')
legend('gen. spinodal law', ['Coulomb' char(39) 's law'])
ylabel('F_f [N]'); xlabel('t [s]');

figure
plot(spinodal.t,spinodal.Fsum, '-',...
     coulomb.t,coulomb.Fsum, '-','LineWidth', 2, 'MarkerSize', 15)
title('net force')
legend('gen. spinodal law', ['Coulomb' char(39) 's law'])
ylabel('F_s [N]'); xlabel('t [s]');
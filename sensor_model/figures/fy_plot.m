clear all; close all

epsilon = 0.01;
x = -0.1:0.001:0.5;
y = epsilon .* exp(-x./epsilon);
figure
plot([0.5 0 0], [0 0 10],'LineWidth', 2)
hold on
plot(x,y,'LineWidth', 2)
axis([-0.1 0.5 0 1])
legend('exact value', 'approximation')
xlabel('$\delta_i$', 'Interpreter', 'latex','FontSize', 16, 'FontWeight', 'bold')
ylabel('$F_Y^i$', 'Interpreter', 'latex','FontSize', 16, 'FontWeight', 'bold')
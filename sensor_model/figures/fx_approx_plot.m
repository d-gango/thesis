clear all; close all

c = 100;
x = -1:0.001:1;
y = tanh(c*x);
figure
plot([-1 0 0 1], [-1 -1 1 1],'LineWidth', 2)
hold on
plot(x,y,'LineWidth', 2)
legend('sign(v)', 'tanh(100v)')
xlabel('v')
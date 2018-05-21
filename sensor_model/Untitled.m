x0 = [0 3 5 7 8 14];
y0 = [10 3 5 -4 7 1];
z0 = rand(1,length(x0))*10;
[x,y] = meshgrid(0:0.1:15);
z = griddata(x0, y0, z0, x, y, 'natural');
% y = interp2(x0, y0, z0, x, x, 'makima');
figure
plot3(x0,y0, z0, '.r', 'MarkerSize', 10)
grid on
hold on
mesh(x,y,z);


function f = lagrange(x0,y0)
syms x
P = 0;
for i = 1:length(x0)
    num = 1;
    den = 1;
    for j = 1:length(x0)
        if j ~= i
           num = num * (x-x0(j));
           den = den * (x0(i) - x0(j));
        else
            num = num * y0(j);
        end
    end
    P = P + num/den;
end
f = matlabFunction(P);
end
load d_F_20.mat
global d contacts n D L h
index = 60;
d = solutions(index).d;
contacts = solutions(index).contacts;
n = params.n;
D = params.D;
h = params.h;
L = D*sin(pi/(2*n));

drawsensor(solutions(index).sol);
title('Solution')
drawsensor(solutions(index).approximation);
title('Approximation');
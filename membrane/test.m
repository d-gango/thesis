load d_F_20.mat
global d contacts n D L h
index = 10;
d = solutions(index).d;
contacts = solutions(index).slip_contacts;
n = params.n;
D = params.D;
h = params.h;
L = D*sin(pi/(2*n));

drawsensor(solutions(index).slip_solution);
title('Solution')
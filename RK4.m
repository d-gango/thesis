function [x] = RK4(fun, x0, t, h)

k1 = fun(x0,          t);
k2 = fun(x0 + h*k1/2, t + h/2);
k3 = fun(x0 + h*k2/2, t + h/2);
k4 = fun(x0 + h*k3  , t + h);

x = x0 + h/6*(k1 + 2*k2 + 2*k3 + k4);



end


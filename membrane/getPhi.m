function phi = getPhi(sol)
par = param();
n = par.n;
phi = zeros(1,n+1);
for i = 1:n+1
    index = (i-1)*3;
    phi(i) = sol(index+1);
end
end
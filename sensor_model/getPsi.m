function psi = getPsi(phi)
n = length(phi);
psi = zeros(1,n);
for j = 1:n
    psi(j) = sum(phi(1:j));
end
end
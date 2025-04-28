function f = GenereteInitialDistribution(n, VT, v_grid, Nv)
f = zeros(Nv,Nv,Nv);
factor = n / (pi^(3/2) * VT^3);
for k = 1:Nv
    for l = 1:Nv
        for m = 1:Nv
            f(m,l,k) = exp(- (v_grid(m)^2 + v_grid(l)^2 + v_grid(k)^2) / VT^2);
        end
    end
end
f = f * factor;
end


function f = GenereteMaxwellDistribution(n, VT, v_grid, V)

arguments
    n
    VT
    v_grid
    V (1,3) double = [0, 0, 0]
end

[Vx, Vy, Vz] = ndgrid(v_grid, v_grid, v_grid);
% Pure maxwell distribution
factor = n / (pi^(3/2) * VT^3);
f = factor * exp(-((Vx - V(1)).^2 + (Vy - V(2)).^2 + (Vz - V(3)).^2) / VT^2);

end


function f = GenereteBiMaxwellDistribution(n, VTpar, VTperp, v_grid)

[Vx, Vy, Vz] = ndgrid(v_grid, v_grid, v_grid);
% Pure maxwell distribution
factor = n / (pi^(3/2) * VTpar * VTperp^2);
f = factor * exp(-Vy.^2/ VTpar^2 -(Vx.^2 + Vz.^2) / VTperp^2);

end


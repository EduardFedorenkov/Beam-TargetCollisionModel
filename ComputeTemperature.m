function T = ComputeTemperature(m, c, vGrid, f)
mult = m / (3 * c^2);
[Vx, Vy, Vz] = ndgrid(vGrid, vGrid, vGrid);
Vsqr = Vx.^2 + Vy.^2 + Vz.^2;
T = sum(Vsqr(:) .* f(:), 'all') / sum(f, 'all') * mult;
end


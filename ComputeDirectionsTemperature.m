function T = ComputeDirectionsTemperature(m, c, vGrid, f, Vmean)
mult = m / (3 * c^2);
[Vx, Vy, Vz] = ndgrid(vGrid, vGrid, vGrid);
Vxsqr = Vx.^2;
Vysqr = Vy.^2;
Vzsqr = Vz.^2;
Tx = (sum(Vxsqr(:) .* f(:), 'all') / sum(f, 'all') - Vmean(1)^2) * mult;
Ty = (sum(Vysqr(:) .* f(:), 'all') / sum(f, 'all') - Vmean(2)^2) * mult;
Tz = (sum(Vzsqr(:) .* f(:), 'all') / sum(f, 'all') - Vmean(3)^2) * mult;
T = [Tx, Ty, Tz];
end


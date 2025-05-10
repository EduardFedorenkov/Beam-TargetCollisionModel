function vel = ComputeVel(vGrid, f)
[Vx, Vy, Vz] = ndgrid(vGrid, vGrid, vGrid);
velx = sum(Vx(:) .* f(:), 'all');
vely = sum(Vy(:) .* f(:), 'all');
velz = sum(Vz(:) .* f(:), 'all');
vel = [velx, vely, velz] / sum(f, 'all');
end


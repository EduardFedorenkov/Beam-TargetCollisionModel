function f = test_nuSourceGenerateDistribution(n, Rp, vGrid)
[Vx, Vy, Vz] = ndgrid(vGrid, vGrid, vGrid);
V = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
f = (V <= Rp) * n / (4 / 3 * pi * Rp^3);
end
function den = ComputeDensity(vGrid, f)
dv = vGrid(2) - vGrid(1);
den = sum(f, 'all') * dv^3;
end


function distr = TargetDistributionFunction(vx, vy, vz, nt, VTt, Vt)
factor = nt / (pi^(3/2) * VTt^3);
distr = factor * exp(- ((vx - Vt(1)).^2 + (vy - Vt(2)).^2 + (vz - Vt(3)).^2) ./ VTt^2);
end


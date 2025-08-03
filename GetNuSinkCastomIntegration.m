function nuSink = GetNuSinkCastomIntegration(Vi_list, Nv, np, VTp, Vp, Rp, crossSection)
Ntotal = size(Vi_list, 1);
nuSink = zeros(Ntotal, 1);
vGridx = linspace(Vp(1) - Rp, Vp(1) + Rp, Nv);
vGridy = linspace(Vp(2) - Rp, Vp(2) + Rp, Nv);
vGridz = linspace(Vp(3) - Rp, Vp(3) + Rp, Nv);
dv = vGridx(2) - vGridx(1);
[Vx, Vy, Vz] = ndgrid(vGridx, vGridy, vGridz);

mask = (Vx - Vp(1)).^2 + (Vy - Vp(2)).^2 + (Vz - Vp(3)).^2 <= Rp^2;

for i = 1:Ntotal
    Vi = Vi_list(i, :);

    integrand = sqrt( (Vi(1) - Vx).^2 + (Vi(2) - Vy).^2 + (Vi(3) - Vz).^2 ) .* ...
        TargetDistributionFunction(Vx, Vy, Vz, np, VTp, Vp);

    integrand(~mask) = 0;

    nuSink(i) = sum(integrand, 'all');
end

nuSink = nuSink * dv^3 * crossSection;
end

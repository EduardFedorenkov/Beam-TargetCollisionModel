function nuSink = GetNuSink(Vi_list, np, VTp, Vp, Rp, crossSection)
Ntotal = size(Vi_list, 1);
nuSink = zeros(Ntotal, 1);
for i = 1:Ntotal
    Vi = Vi_list(i, :);
    
    vxMin = Vp(1) - Rp;
    vyMin = Vp(2) - Rp;
    vzMin = Vp(3) - Rp;

    vxMax = Vp(1) + Rp;
    vyMax = Vp(2) + Rp;
    vzMax = Vp(3) + Rp;

    integrand = @(vx, vy, vz) sqrt( (Vi(1) - vx).^2 + (Vi(2) - vy).^2 + (Vi(3) - vz).^2 ) .* ...
        TargetDistributionFunction(vx, vy, vz, np, VTp, Vp);
    nuSink(i) = crossSection * integral3(integrand, vxMin, vxMax, vyMin, vyMax, vzMin, vzMax);
end
end


function nuSink = GetNuSink(Nv, vGrid, np, VTp, Vp, crossSection)
nuSink = zeros(Nv*Nv*Nv, 1);
for i = 1:Nv*Nv*Nv
    [ki, li, mi] = ind2sub([Nv, Nv, Nv], i);
    Vi = [vGrid(ki), vGrid(li), vGrid(mi)];
    
    vxMin = Vp(1) - 4 / sqrt(2) * VTp;
    vyMin = Vp(2) - 4 / sqrt(2) * VTp;
    vzMin = Vp(3) - 4 / sqrt(2) * VTp;

    vxMax = Vp(1) + 4 / sqrt(2) * VTp;
    vyMax = Vp(2) + 4 / sqrt(2) * VTp;
    vzMax = Vp(3) + 4 / sqrt(2) * VTp;

    integrand = @(vx, vy, vz) sqrt( (Vi(1) - vx).^2 + (Vi(2) - vy).^2 + (Vi(3) - vz).^2 ) .* ...
        TargetDistributionFunction(vx, vy, vz, np, VTp, Vp);
    nuSink(i) = crossSection * integral3(integrand, vxMin, vxMax, vyMin, vyMax, vzMin, vzMax);
end
end


function f = IntegrandFunction(vx, vy, vz, normUji, Vj, mg, mp, np, VTp, Vp)
normU = sqrt( (Vj(1) - vx).^2 + (Vj(2) - vy).^2 + (Vj(3) - vz).^2 );
M = mg + mp;
chi = acos(1 - 2 * (M / (2 * mp) * normUji ./ normU).^2);
if (any(any(1 - 2 * (M / (2 * mp) * normUji ./ normU).^2 > 1) == true))
    disp("cos > 1")
    disp(normU)
    disp(normUji)
end
if (any(any(1 - 2 * (M / (2 * mp) * normUji ./ normU).^2 < -1) == true))
    disp("cos < -1")
    disp(normU)
    disp(normUji)
end
if (any(any(normUji > normU) == true))
    disp("Uji > U")
end
diffCrossSection = GetDiffCrossSection(mg, mp, normU, chi, 1, 1);
f = TargetDistributionFunction(vx, vy, vz, np, VTp, Vp) .* diffCrossSection;
end

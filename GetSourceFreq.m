function sourceFreq = GetSourceFreq(mg, mp, Vi, Vj, np, VTp, Vp, diffCrossSection)
% for now Rp = 4VTp. It suppose that fp - is the Maxwellian distribution
% function. 4 / sqrt(2) > 3 / sqrt(2) - 3 sigma Gause distrubution
[isIntersect, circle] = SourceIntegrationCircle(mg, mp, Vi, Vj, Vp, 4 / sqrt(2) * VTp);
if (isIntersect)
    M = mg + mp;

    integrand = @(t, r) IntegrandSource( ...
    circle.center(1) + r .* cos(t) * circle.ex(1) + r .* sin(t) * circle.ey(1), ...
    circle.center(2) + r .* cos(t) * circle.ex(2) + r .* sin(t) * circle.ey(2), ...
    circle.center(3) + r .* cos(t) * circle.ex(3) + r .* sin(t) * circle.ey(3), ...
    M, mp, Vi, Vj, np, VTp, Vp) .* r;
    sourceFreq = M / (2 * mp) * diffCrossSection * integral2(integrand, 0, 2*pi, 0, circle.r);
else
    sourceFreq = 0;
end
end

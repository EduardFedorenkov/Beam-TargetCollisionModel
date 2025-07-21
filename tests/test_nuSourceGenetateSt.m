function sourceFreq = test_nuSourceGenetateSt(mg, mp, Vi, Vj, np, Rp)
% for now Rp = 4VTp. It suppose that fp - is the Maxwellian distribution
% Rp > 3 / sqrt(2) - 3 sigma Gause distrubution
[isIntersect, circle] = SourceIntegrationCircle(mg, mp, Vi, Vj, [0, 0, 0], Rp);
if (isIntersect)
    integrand = @(t, r) double((circle.center(1) + r .* cos(t) * circle.ex(1) + r .* sin(t) * circle.ey(1)).^2 ...
                             + (circle.center(2) + r .* cos(t) * circle.ex(2) + r .* sin(t) * circle.ey(2)).^2 ...
                             + (circle.center(3) + r .* cos(t) * circle.ex(3) + r .* sin(t) * circle.ey(3)).^2 <= Rp^2) ...
                             .* r .* np / (4/3 * pi * Rp^3);
    sourceFreq = integral2(integrand, 0, 2*pi, 0, circle.r);
else
    sourceFreq = 0;
end
end


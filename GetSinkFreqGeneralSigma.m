function sinkFreq = GetSinkFreqGeneralSigma(mg, mp, Vi, Vj, normUji, np, VTp, Vp, Rp)
% Rp > 3 / sqrt(2) - 3 sigma Gause distrubution
[isIntersect, circle] = SinkIntegrationCircle(mg, mp, Vi, Vj, Vp, Rp);
if (isIntersect)
    integrand = @(t, r) IntegrandFunction( ...
    circle.center(1) + r .* cos(t) * circle.ex(1) + r .* sin(t) * circle.ey(1), ...
    circle.center(2) + r .* cos(t) * circle.ex(2) + r .* sin(t) * circle.ey(2), ...
    circle.center(3) + r .* cos(t) * circle.ex(3) + r .* sin(t) * circle.ey(3), ...
    normUji, Vi, mg, mp, np, VTp, Vp) .* r;
    sinkFreq = integral2(integrand, 0, 2*pi, 0, circle.r);
else
    sinkFreq = 0;
end
end

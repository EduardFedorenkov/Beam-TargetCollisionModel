function f = IntegrandSource(Vpx, Vpy, Vpz, M, mp, normUji, Vj, np, VTp, Vp)
% compute normU = |Vj - Vp|
normU = sqrt( (Vj(1) - Vpx).^2 + (Vj(2) - Vpy).^2 + (Vj(3) - Vpz).^2 );

% det(J) Jacobian
J = 2 * M / mp * normU.^2 / normUji^2;
f = normUji ./ normU.^2 .* J .* TargetDistributionFunction(Vpx, Vpy, Vpz, np, VTp, Vp);
end

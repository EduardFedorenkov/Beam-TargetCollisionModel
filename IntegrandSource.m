function f = IntegrandSource(Vpx, Vpy, Vpz, M, mp, Vi, Vj, np, VTp, Vp)
% compute normU = |Vj - Vp|
normU = sqrt( (Vj(1) - Vpx).^2 + (Vj(2) - Vpy).^2 + (Vj(3) - Vpz).^2 );

% compute Uji
Uji = Vj - Vi;
normUji = norm(Uji);

% det(J) Jacobian
J = 1 - M/mp + 2 * mp / M * normU.^2 / normUji^2;
f = J .* normU .* TargetDistributionFunction(Vpx, Vpy, Vpz, np, VTp, Vp);
end

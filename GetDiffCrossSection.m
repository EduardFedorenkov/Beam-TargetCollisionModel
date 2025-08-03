function diffSigma = GetDiffCrossSection(mg, mp, u, chi, Z1, Z2)
electronCharge = 4.80320425e-10;                    % СГС
eVtoErg = 1.602176634e-12;                          % eV to Erg convertion
c = 2.99792458e10;                                  % Speed of light [cm/s]
mu = mg * mp / (mg + mp) * eVtoErg / c^2;
diffSigma = Z1 * Z2 * electronCharge^4 ./ (16 * mu^2 .* u.^4) ./ sin(chi / 2).^4;
end


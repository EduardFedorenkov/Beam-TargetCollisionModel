function nu = PlasmaIonsTermolization(ma, mb, Za, Zb, na, nb, Ta, Tb)
%% Set constants
% Physical consts
factor = 1.8e-19;
c = 2.99792458e10;                                  % Speed of light [cm/s]
eVtoErg = 1.602176634e-12;                          % Convertion coef from [eV] to [Erg]
mp = 938.272e6;
mua = ma / mp;
mub = mb / mp;
maa = ma * eVtoErg / c^2;
mbb = mb * eVtoErg / c^2;
lambda = 23 - log(Za * Zb * (mua + mub) / (mua * Tb + mub * Ta) * sqrt(na * Za^2 / Ta + nb * Zb^2 / Tb));
% lambda = 10;

nu = factor * sqrt(maa * mbb) * Za^2 * Zb^2 * nb * lambda / (maa * Tb + mbb * Ta)^(3/2);
end


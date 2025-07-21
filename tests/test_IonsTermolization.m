%% Set constants
% Physical consts
c = 2.99792458e10;                                  % Speed of light [cm/s]
eVtoErg = 1.602176634e-12;                          % Convertion coef from [eV] to [Erg]
m = 938.272e6;                                      % Mass of proton [eV]

% B ion
Zb = 1;
mb = m;                                             % Ions mass [eV]
nb = 1e14;                                          % Ions density [cm^{-3}]
Tb = 5;                                             % Ions temperature [eV]
VTb = sqrt(2 * Tb / mb) * c;                        % Ions termal vel [cm /s]
Vb = [0, 0, 0];                                     % Ions vel [cm / s] (vectro size of 3!!!)
Rb = 5 / sqrt(2) * VTb;                             % 5-siga maxwell range

% A ion
Za = 1;
ma = m;                                             % Gas mass [eV]
na = 1e14;                                          % Gas dencity [cm^{-3}]
Ta = 3;                                             % Gas temperature [eV]
VTa = sqrt(2 * Ta / ma) * c;                        % Gas termal vel [cm /s]

timeEstimation = 1 / PlasmaIonsTermolization(ma, mb, Za, Zb, na, nb, Ta, Tb);
dt = timeEstimation / 20;
N = 1000;

TaVec = TermalizationSolver(ma, mb, Za, Zb, na, nb, Ta, Tb, dt, N);

time = 0:dt:(N * dt);
plot(time, TaVec);
title('T(t) reference');
xlabel('time [s]');
ylabel('T [eV]');

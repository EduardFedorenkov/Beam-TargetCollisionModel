%% Set constants
% Physical consts
c = 3 * 10^10;                                      % Speed of light [cm/s]
eVtoErg = 1.6e-12;                                  % Convertion coef from [eV] to [Erg]
m = 938.27 * 10^6;                                  % Mass of proton [eV]
eps = 3 / sqrt(2);                                  % 3-siga maxwell range

% Plasma parameters
mp = 2 * m;                                         % Ions mass [eV]
np = 1e14;                                          % Ions density [cm^{-3}]
Tp = 50;                                            % Ions temperature [eV]
VTp = sqrt(2 * Tp / mp) * c;                        % Ions termal vel [cm /s]
Vp = [0, sqrt(2 * 200 / mp) * c, 0];                % Ions vel [cm / s] (vectro size of 3!!!)

% Gas parameters
mg = 4 * m;                                         % Gas mass [eV]
ng = 1e14;                                          % Gas dencity [cm^{-3}]
Tg = 1;                                             % Gas temperature [eV]
VTg = sqrt(2 * Tg / mg) * c;                        % Gas termal vel [cm /s]

% GAS-ION Cross section
diffCrossSection = 1e-16 / (4 * pi);                % HH collisions [cm^2]
crossSection = 1e-16;

%% Set model parameters
Nv = 11;
vGrid = linspace(-2 * eps * VTg, 2 * eps * VTg, Nv);
dv = vGrid(2) - vGrid(1);

%% Begin computation
fg = GenereteInitialDistribution(ng, VTg, vGrid, Nv);

nuSink = GetNuSink(Nv, vGrid, np, VTp, Vp, crossSection);
nuSource = GetNuSource(Nv, vGrid, mg, mp, np, VTp, Vp, diffCrossSection) * dv^3;

st = GetSt(nuSink, nuSource, fg, Nv);

%% Plot results
pcolor(vGrid, vGrid, st(:,:,6)');
shading flat;
shading interp;
title('St, v_z = 0');
xlabel('v_x [cm / s]'); 
ylabel('v_y [cm / s]');

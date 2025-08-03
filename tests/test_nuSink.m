%% Set constants
% Physical consts
c = 2.99792458e10;                                  % Speed of light [cm/s]
eVtoErg = 1.602176634e-12;                          % Convertion coef from [eV] to [Erg]
m = 938.272e6;                                      % Mass of proton [eV]

% Plasma parameters
mp = m;                                             % Ions mass [eV]
np = 1e14;                                          % Ions density [cm^{-3}]
Tp = 5;                                             % Ions temperature [eV]
VTp = sqrt(2 * Tp / mp) * c;                        % Ions termal vel [cm /s]
Vp = [0, 0, 0];                                     % Ions vel [cm / s] (vectro size of 3!!!)
Rp = 5 / sqrt(2) * VTp;                             % 5-siga maxwell range

% Gas parameters
mg = m;                                             % Gas mass [eV]
ng = 1e14;                                          % Gas dencity [cm^{-3}]
Tg = 5;                                             % Gas temperature [eV]
VTg = sqrt(2 * Tg / mg) * c;                        % Gas termal vel [cm /s]

% GAS-ION Cross section
diffCrossSection = 1e-16 / (4 * pi);                % HH collisions [cm^2]
crossSection = 1e-16;

%% Set model parameters
Nv = 11;
Ntotal = Nv^3;
vGrid = linspace(-Rp, Rp, Nv);
dv = vGrid(2) - vGrid(1);

[Vx, Vy, Vz] = ndgrid(vGrid, vGrid, vGrid);

%% Begin computation
fg = GenereteInitialDistribution(ng, VTg, vGrid, Nv);

% Precompute all velocity vectors
Vi_list = zeros(Ntotal, 3);
for i = 1:Ntotal
    [ki, li, mi] = ind2sub([Nv, Nv, Nv], i);
    Vi_list(i, :) = [vGrid(ki), vGrid(li), vGrid(mi)];
end

% nuSink = GetNuSink(Vi_list, np, VTp, Vp, Rp, crossSection);
nuSink = GetNuSinkCastomIntegration(Vi_list, Nv, np, VTp, Vp, Rp, crossSection);

stSink = nuSink;
stSink = reshape(stSink, [Nv, Nv, Nv]);

V = sqrt(Vx.^2 + Vy.^2 + Vz.^2) / VTp;

stSinkTest = crossSection * np * VTp * ( ( 1 + 2 * V.^2 ) ./ (2 * V ) .* erf(V) + 1/sqrt(pi) * exp(-V.^2) );
stSinkTest(isnan(stSinkTest)) = 2 / sqrt(pi) * crossSection * np * VTp; % in zero point

test = (stSink - stSinkTest) ./ stSinkTest * 100;

figure(1);
pcolor(vGrid, vGrid, test(:,:,6)');
shading flat;
shading interp;
title('(stSink - stSinkTest) / stSinkTest * 100 [%]');
xlabel('v_x [cm / s]'); 
ylabel('v_y [cm / s]');


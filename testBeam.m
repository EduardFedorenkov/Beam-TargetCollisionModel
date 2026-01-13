%% Set constants
% Physical consts
c = 2.99792458e10;                                  % Speed of light [cm/s]
eVtoErg = 1.602176634e-12;                          % Convertion coef from [eV] to [Erg]
m = 938.272e6;                                      % Mass of proton [eV]

% Plasma parameters
mp = m;                                             % Ions mass [eV]
np = 5.5e13;                                        % Ions density [cm^{-3}]
Tp = 2300;                                          % Ions temperature [eV]
VTp = sqrt(2 * Tp / mp) * c;                        % Ions termal vel [cm /s]
Ep = 60000;                                         % Ion beam kinetic enegry [eV]
Vp = [0, -sqrt(2 * Ep / mp) * c, 0];                % Ions vel [cm / s] (vectro size of 3!!!)
Rp = 5 / sqrt(2) * VTp;                             % Rp-siga maxwell range

% Gas parameters
mg = m;                                             % Gas mass [eV]
ng = 1e13;                                          % Gas dencity [cm^{-3}]
TgPar = Ep * 1e-3;                                  % Gas parallel temperature [eV]
TgPerp = 30;                                        % Gas perpendicular temperature [eV]
VTgPar = sqrt(2 * TgPar / mg) * c;                  % Gas parallel termal vel [cm /s]
VTgPerp = sqrt(2 * TgPerp / mg) * c;                % Gas perpendicular termal vel [cm /s]

% GAS-ION Cross section
M = mg + mp;

%% Set model parameters
Nv = 13;
Ntotal = Nv^3;
Vmax = 5 / sqrt(2) * sqrt(2 * 200 / mg) * c;
vGrid = linspace(-Vmax, Vmax, Nv);
dv = vGrid(2) - vGrid(1);

% Precompute all velocity vectors
Vi_list = zeros(Ntotal, 3);
for i = 1:Ntotal
    [ki, li, mi] = ind2sub([Nv, Nv, Nv], i);
    Vi_list(i, :) = [vGrid(ki), vGrid(li), vGrid(mi)];
end

% Precompute normUji matrix using vectorization
Vi_3d = reshape(Vi_list, [Ntotal, 1, 3]);
Vj_3d = reshape(Vi_list, [1, Ntotal, 3]);
diff = Vj_3d - Vi_3d;
normUji = sqrt(sum(diff.^2, 3));

% Compute 1 ./ normUji with zero diagonal elements
normUjiInv = 1 ./ normUji;
normUjiInv(logical(eye(size(normUji)))) = zeros(Ntotal, 1);

%% Begin computation
% Test for general cross section:
sigmaSourceFactor = (M / mp)^2 * dv^3 .* normUjiInv;
nuSource = sigmaSourceFactor .* GetNuSourceGeneralSigma(Vi_list, normUji, mg, mp, np, VTp, Vp, Rp);
nuSinkMat = sigmaSourceFactor .* GetNuSinkMatGeneralSigma(Vi_list, normUji, mg, mp, np, VTp, Vp, Rp);

nuSink = sum(nuSinkMat, 2);

fg = GenereteBiMaxwellDistribution(ng, VTgPar, VTgPerp, vGrid);
st = GetSt(nuSink, nuSource, fg, Nv);

%% Evolution of fg over time
Nl = 5000;
Vmean = sqrt(2 * Ep / mp) * c;
dl = 0.01;

n = zeros(Nl+1, 1);
V = zeros(Nl+1, 3);
T = zeros(Nl+1, 3);

n(1) = ComputeDensity(vGrid, fg);
V(1,:) = ComputeVel(vGrid, fg);
T(1,:) = ComputeDirectionsTemperature(mg, c, vGrid, fg, [0, 0, 0]);

f = fg;
for l = 1:Nl
    f = SemiImplicitTransportScheme(f, nuSource, nuSink, dl, Vmean);
    
    n(l + 1) = ComputeDensity(vGrid, f);
    V(l + 1, :) = ComputeVel(vGrid, f);
    T(l + 1, :) = ComputeDirectionsTemperature(mg, c, vGrid, f, V(l + 1, :));
end
f = reshape(f, [Nv, Nv, Nv]);

%% Plot results
figure;
pcolor(vGrid, vGrid, st(:,:,7)');
colormap("jet");
colorbar;
shading flat;
shading interp;
title('St, v_z = 0');
xlabel('v_x [cm / s]'); 
ylabel('v_y [cm / s]');

figure;
pcolor(vGrid, vGrid, fg(:,:,7)');
colormap("jet");
colorbar;
shading flat;
shading interp;
title('Initial distribution, v_z = 0');
xlabel('v_x [cm / s]'); 
ylabel('v_y [cm / s]');

strLenEnd = num2str(dl * Nl);
strNl = num2str(Nl);
figure;
pcolor(vGrid, vGrid, f(:,:,7)');
colormap("jet");
colorbar;
shading flat;
shading interp;
title(['f(t = ', strLenEnd, ' sec ', strNl, ' time steps), v_z = 0']);
xlabel('v_x [cm / s]');
ylabel('v_y [cm / s]');

% Plot moments of gas distribution function
figure;
time = linspace(0, Nl*dl, Nl + 1);
plot(time, n, 'r-', 'LineWidth', 2);
grid on;
title('Gas Density');
xlabel('length [cm]');
ylabel('n [cm^{-3}]');

figure;
plot(time, T(:, 2), 'r-', 'LineWidth', 2);
grid on;
title('Gas Temperature ||');
xlabel('length [cm]');
ylabel('Ty [eV]');

figure;
plot(time, T(:, 1), 'r-', 'LineWidth', 2);
grid on;
title('Gas Temperature \perp');
xlabel('length [cm]');
ylabel('Tx [eV]');

figure;
plot(time, T(:, 3), 'r-', 'LineWidth', 2);
grid on;
title('Gas Temperature \perp');
xlabel('length [cm]');
ylabel('Tz [eV]');

figure;
plot(time, V(:, 2), 'r-', 'LineWidth', 2);
grid on;
title('Gas v_{||}');
xlabel('length [cm]');
ylabel('V [cm / s]');

figure;
plot(time, V(:, 1), 'r-', 'LineWidth', 2);
grid on;
title('Gas v_{\perp}');
xlabel('length [cm]');
ylabel('V [cm / s]');

figure;
plot(time, V(:, 3), 'r-', 'LineWidth', 2);
grid on;
title('Gas v_{\perp}');
xlabel('length [cm]');
ylabel('V [cm / s]');

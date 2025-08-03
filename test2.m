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
Vp = [0, VTp, 0];                                   % Ions vel [cm / s] (vectro size of 3!!!)
Rp = 5 / sqrt(2) * VTp;                             % Rp-siga maxwell range

% Gas parameters
mg = m;                                             % Gas mass [eV]
ng = 1e14;                                          % Gas dencity [cm^{-3}]
Tg = 3;                                             % Gas temperature [eV]
VTg = sqrt(2 * Tg / mg) * c;                        % Gas termal vel [cm /s]

% GAS-ION Cross section
M = mg + mp;
diffCrossSection = 1e-16 / (4 * pi);                % HH collisions [cm^2]
crossSection = 1e-16;

%% Set model parameters
Nv = 13;
Ntotal = Nv^3;
Vmax = 5 / sqrt(2) * sqrt(2 * Tp / mg) * c;
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
% Test with constant cross section:
constSigmaSourceFactor = (M / mp)^2 * diffCrossSection * dv^3 .* normUjiInv;
nuSource = constSigmaSourceFactor .* GetNuSourceConstSigma(Vi_list, mg, mp, np, VTp, Vp, Rp);
nuSinkMat = constSigmaSourceFactor .* GetNuSinkMatConstSigma(Vi_list, mg, mp, np, VTp, Vp, Rp);
nuSinkOld = GetNuSink(Vi_list, np, VTp, Vp, Rp, crossSection);

% Test for general cross section:
% constSigmaSourceFactor = (M / mp)^2 * dv^3 .* normUjiInv;
% nuSource = constSigmaSourceFactor .* GetNuSourceGeneralSigma(Vi_list, normUji, mg, mp, np, VTp, Vp, Rp);
% nuSinkMat = constSigmaSourceFactor .* GetNuSinkMatGeneralSigma(Vi_list, normUji, mg, mp, np, VTp, Vp, Rp);

nuSink = sum(nuSinkMat, 2);
%nuSink = nuSinkOld;

fg = GenereteMaxwellDistribution(ng, VTg, vGrid);
fgFinal = GenereteMaxwellDistribution(ng, VTp * sqrt(mp / mg), vGrid, Vp * sqrt(mp / mg));

st = GetSt(nuSink, nuSource, fg, Nv);

figure(1);
pcolor(vGrid, vGrid, st(:,:,6)');
colormap("jet");
colorbar;
shading flat;
shading interp;
title('St, v_z = 0');
xlabel('v_x [cm / s]'); 
ylabel('v_y [cm / s]');

%% Evolution of fg over time
Nt = 400;
dt = 1e-6;

n = zeros(Nt+1, 1);
V = zeros(Nt+1, 3);
T = zeros(Nt+1, 1);
meanDistrDiviation = zeros(Nt+1, 1);

n(1) = ComputeDensity(vGrid, fg);
V(1,:) = ComputeVel(vGrid, fg);
T(1) = ComputeTemperature(mg, c, vGrid, fg, [0, 0, 0]);
meanDistrDiviation(1) = mean(abs(fg(:) - fgFinal(:)) ./ fgFinal(:)) * 100;

f = fg;
for t = 1:Nt
    f = SemiImplicitTimeScheme(f, nuSource, nuSink, dt);
    meanDistrDiviation(t + 1) = mean(abs(f(:) - fgFinal(:)) ./ fgFinal(:)) * 100;
    
    n(t + 1) = ComputeDensity(vGrid, f);
    V(t + 1, :) = ComputeVel(vGrid, f);
    T(t + 1) = ComputeTemperature(mg, c, vGrid, f, V(t + 1, :));
end
f = reshape(f, [Nv, Nv, Nv]);

%% Plot results
figure(2);
pcolor(vGrid, vGrid, fg(:,:,7)');
colormap("jet");
colorbar;
shading flat;
shading interp;
title('Initial distribution, v_z = 0');
xlabel('v_x [cm / s]'); 
ylabel('v_y [cm / s]');

strTimeEnd = num2str(dt * Nt);
strNt = num2str(Nt);
figure(3);
pcolor(vGrid, vGrid, f(:,:,7)');
colormap("jet");
colorbar;
shading flat;
shading interp;
title(['f(t = ', strTimeEnd, ' sec ', strNt, ' time steps), v_z = 0']);
xlabel('v_x [cm / s]');
ylabel('v_y [cm / s]');

%% Plot moments of gas distribution function
figure(4)
time = linspace(0, Nt*dt, Nt + 1);
plot(time, nnew, 'r-', 'LineWidth', 2);
hold on;
plot(time, n, 'b--', 'LineWidth', 2);
grid on;
title('Gas Density');
xlabel('time [s]');
ylabel('n [cm^{-3}]');

figure(5)
plot(time, Tnew, 'r-', 'LineWidth', 2);
hold on;
plot(time, T, 'b--', 'LineWidth', 2);
grid on;
title('Gas Temperature');
xlabel('time [s]');
ylabel('T [eV]');

figure(6)
plot(time, Vnew(:, 2), 'r-', 'LineWidth', 2);
hold on;
plot(time, V(:, 2), 'b--', 'LineWidth', 2);
grid on;
title('Gas v_{||}');
xlabel('time [s]');
ylabel('V [cm / s]');

figure(7)
plot(time, V(:, 1), 'r-', 'LineWidth', 2);
grid on;
title('Gas v_{\perp}');
xlabel('time [s]');
ylabel('V [cm / s]');

figure(8)
plot(time, V(:, 3), 'r-', 'LineWidth', 2);
grid on;
title('Gas v_{\perp}');
xlabel('time [s]');
ylabel('V [cm / s]');

figure(9)
plot(time, meanDistrDiviation, 'g-', 'LineWidth', 2);
grid on;
title('mean deviation from the stationary distribution: |f - f_{exact}| / f_{exact} * 100 %');
xlabel('time [s]');
ylabel('<\delta f / f> [%]');

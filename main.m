%% Set constants
% Physical consts
c = 2.99792458e10;                                  % Speed of light [cm/s]
eVtoErg = 1.602176634e-12;                          % Convertion coef from [eV] to [Erg]
m = 938.272e6;                                      % Mass of proton [eV]
eps = 3 / sqrt(2);                                  % 3-siga maxwell range

% Plasma parameters
mp = m;                                             % Ions mass [eV]
np = 1e14;                                          % Ions density [cm^{-3}]
Tp = 3;                                             % Ions temperature [eV]
VTp = sqrt(2 * Tp / mp) * c;                        % Ions termal vel [cm /s]
Vp = [0, 0, 0];                                     % Ions vel [cm / s] (vectro size of 3!!!)

% Gas parameters
mg = m;                                             % Gas mass [eV]
ng = 1e14;                                          % Gas dencity [cm^{-3}]
Tg = 1;                                             % Gas temperature [eV]
VTg = sqrt(2 * Tg / mg) * c;                        % Gas termal vel [cm /s]

% GAS-ION Cross section
diffCrossSection = 1e-16 / (4 * pi);                % HH collisions [cm^2]
crossSection = 1e-16;

%% Set model parameters
Nv = 11;
Ntotal = Nv^3;
vGrid = linspace(-eps * VTp, eps * VTp, Nv);
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
normUji_matrix = sqrt(sum(diff.^2, 3));

%% Begin computation
nuSource = GetNuSource(Vi_list, normUji_matrix, mg, mp, np, VTp, Vp, diffCrossSection) * dv^3;
nuSink = GetNuSink(Vi_list, np, VTp, Vp, crossSection);

fg = GenereteInitialDistribution(ng, VTg, vGrid, Nv);

st = GetSt(nuSink, nuSource, fg, Nv);

figure(1);
pcolor(vGrid, vGrid, st(:,:,6)');
shading flat;
shading interp;
title('St, v_z = 0');
xlabel('v_x [cm / s]'); 
ylabel('v_y [cm / s]');

%% Evolution of fg over time
Nt = 300;
dt = 1e-6;

n = zeros(Nt+1, 1);
V = zeros(Nt+1, 3);
T = zeros(Nt+1, 1);

n(1) = ComputeDensity(vGrid, fg);
V(1,:) = ComputeVel(vGrid, fg);
T(1) = ComputeTemperature(mg, c, vGrid, fg);

f = fg;
for t = 1:Nt
    f = SemiImplicitTimeScheme(f, nuSource, nuSink, dt);

    n(t + 1) = ComputeDensity(vGrid, f);
    V(t + 1, :) = ComputeVel(vGrid, f);
    T(t + 1) = ComputeTemperature(mg, c, vGrid, f);
end
f = reshape(f, [Nv, Nv, Nv]);

%% Plot results
figure(2);
pcolor(vGrid, vGrid, fg(:,:,6)');
shading flat;
shading interp;
title('Initial distribution, v_z = 0');
xlabel('v_x [cm / s]'); 
ylabel('v_y [cm / s]');

strTimeEnd = num2str(dt * Nt);
strNt = num2str(Nt);
figure(3);
pcolor(vGrid, vGrid, f(:,:,6)');
shading flat;
shading interp;
title(['f(t = ', strTimeEnd, ' sec ', strNt, ' time steps), v_z = 0']);
xlabel('v_x [cm / s]');
ylabel('v_y [cm / s]');

%% Plot moments of gas distribution function
figure(4)
time = linspace(0, Nt*dt, Nt + 1);
plot(time, n);
title('Gas Density');
xlabel('time [s]');
ylabel('n [cm^{-3}]');

figure(5)
plot(time, T);
title('Gas Temperature');
xlabel('time [s]');
ylabel('T [eV]');

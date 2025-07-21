%% Set constants
% Physical consts
c = 2.99792458e10;                                  % Speed of light [cm/s]
eVtoErg = 1.602176634e-12;                          % Convertion coef from [eV] to [Erg]
m = 938.272e6;                                      % Mass of proton [eV]

% Gas and Plasma parameters
mp = m;                                             % Ions mass [eV]
mg = m;                                             % Gas mass [eV]
n = 1e14;                                           % density [cm^{-3}]
T = 1;                                              % Ions temperature [eV]
VT = sqrt(2 * T / m) * c;                           % Ions termal vel [cm /s]

% Test scenario: fg = fp = n / (4/3piR^3) = const inside the ball B(0, R)
R = 3 / sqrt(2) * VT;

M = mg + mp;

% GAS-ION Cross section
crossSection = 1e-16;
diffCrossSection = crossSection / (4 * pi);                % HH collisions [cm^2]


%% Set model parameters
Nv = 11;
Ntotal = Nv^3;
vGrid = linspace(-1.5 * R, 1.5 * R, Nv);
dv = vGrid(2) - vGrid(1);

%% Work with distribution function
fg = test_nuSourceGenerateDistribution(n, R, vGrid);
fp = fg;

% Set Vi, Vj and Uji
i = [4, 5, 6];
Vi = [vGrid(i(1)), vGrid(i(2)), vGrid(i(3))];

stRef = 0;
nu = zeros(Ntotal, 1);
for j = 1:Ntotal
    [jx, jy, jz] = ind2sub([Nv, Nv, Nv], j);
    Vj = [vGrid(jx), vGrid(jy), vGrid(jz)];
    Uji = Vj - Vi;
    normUji = norm(Uji);
    D = ( M / (2 * mp) * normUji^2 - dot(Uji, Vj) ) / normUji;
    if (R >= D) && (norm(Vj) <= R) && (normUji > 0)
        stRef = stRef + 1 / normUji * (R^2 - D^2);
    end
    if (normUji > 0)
        nu(j) = test_nuSourceGenetateSt(mg, mp, Vi, Vj, n, R) / normUji;
    end
end

stRef = stRef * 3^2 / (pi * 4^2) * (M / mp)^2 * diffCrossSection * n^2 / R^6 * dv^3;
nu = nu * (M / mp)^2 * diffCrossSection * dv^3;
st = dot(nu, fg(:));

format longE
disp("source reference : ");
disp(stRef);
disp("source model : ");
disp(st);


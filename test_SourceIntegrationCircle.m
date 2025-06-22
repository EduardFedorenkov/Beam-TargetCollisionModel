%% Test for function SourceIntegrationCircle.m
%% test 1 parameters
mg = 1;
mp = 1;
Vi = [2, 3, 4];
Vj = [-1, -2, -3];
Vp = [-1, 1, 3];
Rp = 10;

% Generate plane
N = Vj - Vi;                                        % Plane normal
D = (mg + mp) / (2 * mp) * norm(N)^2 - dot(N, Vj);  % Plane offset
[x, y] = meshgrid(-2 * Rp:0.1 * 2 * Rp:2 * Rp);     % Generate x and y data
z = -1 / N(3) * (N(1) * x + N(2) * y + D);          % Plane equation

% Generate Sphere
[X, Y, Z] = sphere;                              % Unit Sphere in [0, 0, 0]
X2 = X * Rp + Vp(1);
Y2 = Y * Rp + Vp(2);
Z2 = Z * Rp + Vp(3);

% Plot the results
figure(1)
surf(X2, Y2, Z2, 'EdgeColor', 'none', 'FaceAlpha', 0.4) %plot sphere
hold on
surf(x, y, z, 'EdgeAlpha', 0.4, 'FaceAlpha', 0.4) %Plot the surface

% Compute intersection circle
[isIntersect, circle] = SourceIntegrationCircle(mg, mp, Vi, Vj, Vp, Rp);

% Check asserts
assert(isIntersect == true, "SourceIntegrationCircle return no" + ...
    "intercection with plane!");
assert(circle.r <= Rp + eps, "circle radius r > Rp!");
assert(dot(circle.ex, circle.ey) < 3 * eps, "ex is not orthogonal to ey!");
assert(dot(circle.ex, circle.ex) < 1 + 3 * eps, "|ex| is not 1!");
assert(dot(circle.ey, circle.ey) < 1 + 3 * eps, "|ex| is not 1!");

% Generate edge of thr circle
t = linspace(0, 2 * pi, 100);
xCircle = circle.center(1) + circle.r * cos(t) * circle.ex(1) + circle.r * sin(t) * circle.ey(1);
yCircle = circle.center(2) + circle.r * cos(t) * circle.ex(2) + circle.r * sin(t) * circle.ey(2);
zCircle = circle.center(3) + circle.r * cos(t) * circle.ex(3) + circle.r * sin(t) * circle.ey(3);
plot3(xCircle, yCircle, zCircle, 'r-', 'LineWidth', 2);

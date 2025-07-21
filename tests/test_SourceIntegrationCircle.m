%% Test for function SourceIntegrationCircle.m
%% test 1 parameters
mg = 1;
mp = 1;
Vi = [2, 3, 4];
Vj = [-1, -2, -3];
Vp = [-1, 1, 3];
Rp = 10;

% Generate plane Vp'
N = Vj - Vi;                                            % Plane normal
D1 = (mg + mp) / (2 * mp) * norm(N)^2 - dot(N, Vj);     % Plane offset
[x1, y1] = meshgrid(-2 * Rp:0.1 * 2 * Rp:2 * Rp);       % Generate x and y data
z1 = -1 / N(3) * (N(1) * x1 + N(2) * y1 + D1);          % Plane equation

% Generate plane Vp'
D2 = -(mg + mp) / (2 * mp) * norm(N)^2 - dot(N, Vi);    % Plane offset
[x2, y2] = meshgrid(-2 * Rp:0.1 * 2 * Rp:2 * Rp);       % Generate x and y data
z2 = -1 / N(3) * (N(1) * x2 + N(2) * y2 + D2);          % Plane equation

% Generate line from Vi to Vj
l = linspace(0, 1, 100);
x3 = Vi(1) + l * (Vj(1) - Vi(1));
y3 = Vi(2) + l * (Vj(2) - Vi(2));
z3 = Vi(3) + l * (Vj(3) - Vi(3));

% Generate Sphere
[X, Y, Z] = sphere;                              % Unit Sphere in [0, 0, 0]
X2 = X * Rp + Vp(1);
Y2 = Y * Rp + Vp(2);
Z2 = Z * Rp + Vp(3);

% Plot the results
figure(1)
surf(X2, Y2, Z2, 'EdgeColor', 'none', 'FaceAlpha', 0.1) %plot sphere
hold on
surf(x1, y1, z1, 'EdgeAlpha', 0, 'FaceAlpha', 0.4) %Plot the surface Vp'
surf(x2, y2, z2, 'EdgeAlpha', 0, 'FaceAlpha', 0.4) %Plot the surface Vp
%plot3(x3, y3, z3, 'b-', 'LineWidth', 2) %Plot the line from Vi to Vj

% Compute intersection circle
[isIntersectSource, circleSource] = SourceIntegrationCircle(mg, mp, Vi, Vj, Vp, Rp);
[isIntersectSink, circleSink] = SinkIntegrationCircle(mg, mp, Vi, Vj, Vp, Rp);

% Check asserts
assert(isIntersectSource == true, "SourceIntegrationCircle return no" + ...
    "intercection with plane!");
assert(circleSource.r <= Rp + eps, "circle radius r > Rp!");
assert(dot(circleSource.ex, circleSource.ey) < 3 * eps, "ex is not orthogonal to ey!");
assert(dot(circleSource.ex, circleSource.ex) < 1 + 3 * eps, "|ex| is not 1!");
assert(dot(circleSource.ey, circleSource.ey) < 1 + 3 * eps, "|ex| is not 1!");

% Generate edge of thr circle
t = linspace(0, 2 * pi, 100);
xCircleSource = circleSource.center(1) + circleSource.r * cos(t) * circleSource.ex(1) + circleSource.r * sin(t) * circleSource.ey(1);
yCircleSource = circleSource.center(2) + circleSource.r * cos(t) * circleSource.ex(2) + circleSource.r * sin(t) * circleSource.ey(2);
zCircleSource = circleSource.center(3) + circleSource.r * cos(t) * circleSource.ex(3) + circleSource.r * sin(t) * circleSource.ey(3);
plot3(xCircleSource, yCircleSource, zCircleSource, 'r-', 'LineWidth', 2);

xCircleSink = circleSink.center(1) + circleSink.r * cos(t) * circleSink.ex(1) + circleSink.r * sin(t) * circleSink.ey(1);
yCircleSink = circleSink.center(2) + circleSink.r * cos(t) * circleSink.ex(2) + circleSink.r * sin(t) * circleSink.ey(2);
zCircleSink = circleSink.center(3) + circleSink.r * cos(t) * circleSink.ex(3) + circleSink.r * sin(t) * circleSink.ey(3);
plot3(xCircleSink, yCircleSink, zCircleSink, 'y-', 'LineWidth', 2);

% Generate line from Vi to Vj
x4 = circleSink.center(1) + l * (circleSource.center(1) - circleSink.center(1));
y4 = circleSink.center(2) + l * (circleSource.center(2) - circleSink.center(2));
z4 = circleSink.center(3) + l * (circleSource.center(3) - circleSink.center(3));
plot3(x4, y4, z4, 'g-', 'LineWidth', 2) %Plot the line from C1 to C2

x5 = circleSink.center(1) + l * circleSink.r * circleSink.ex(1);
y5 = circleSink.center(2) + l * circleSink.r * circleSink.ex(2);
z5 = circleSink.center(3) + l * circleSink.r * circleSink.ex(3);
plot3(x5, y5, z5, '-b', 'LineWidth', 2) %Plot the line from C1 to C2

x5 = circleSink.center(1) + l * circleSink.r * circleSink.ey(1);
y5 = circleSink.center(2) + l * circleSink.r * circleSink.ey(2);
z5 = circleSink.center(3) + l * circleSink.r * circleSink.ey(3);
plot3(x5, y5, z5, '-b', 'LineWidth', 2) %Plot the line from C1 to C2

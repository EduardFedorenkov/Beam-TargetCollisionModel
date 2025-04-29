%% Test for function SourceIntegrationCircle.m
% test 1
mg = 1;
mp = 1;
Vi = [2, 3, 4];
Vj = [-1, -2, -3];
Vp = [-1, 1, 3];
Rp = 10;

N = Vj - Vi;
D = (mg + mp) / (2 *mp) * norm(N)^2 - dot(N, Vj);
[x y] = meshgrid(-2 * Rp:0.1:2 * Rp); % Generate x and y data
z = -1/N(3)*(N(1)*x + N(2)*y + D); % Solve for z data

[X,Y,Z] = sphere;
X2 = X * Rp + Vp(1);
Y2 = Y * Rp + Vp(2);
Z2 = Z * Rp + Vp(3);

surf(X2,Y2,Z2) %plot sphere
hold on
surf(x,y,z) %Plot the surface

circle = struct();
[isIntersect, circle] = SourceIntegrationCircle(mg, mp, Vi, Vj, Vp, Rp);
assert(isIntersect == true, "SourceIntegrationCircle return no intercection with plane!");
assert(circle.r <= Rp + eps, "circle radius r > Rp!");
assert(dot(circle.ex, circle.ey) < 3 * eps, "ex is not orthogonal to ey!");

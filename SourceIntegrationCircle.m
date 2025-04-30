function [isIntersect, circle] = SourceIntegrationCircle(mg, mp, Vi, Vj, Vp, Rp)
% we use plane equation: (n * r) + d = 0, where n and r is 3D vectors
% compute plane parameters "n" and "d"
n = Vj - Vi;
normN = norm(n);
d = (mg + mp) / (2 * mg) * normN^2 - dot(n, Vj);

% compute distance between Vp and the plane
dist = abs(dot(n, Vp) + d) / normN;
if (dist > Rp)
    % Sphere don't intersect the plane
    isIntersect = false;
else 
    isIntersect = true;
    circle.center = Vp - (dot(n, Vp) + d) / normN^2 * n;
    circle.r = sqrt(Rp^2 - dist^2);
    
    % compute basis on plane
    % choose any vector \vec{a} not parallel to \vec{n}
    aTmp = eye(3);
    a = [];
    for i = 1:3
        if (norm(cross(aTmp(i,:), n)) > 100 * eps)
            a = aTmp(i,:);
            break
        end
    end
    % compute ex and ey
    circle.ex = cross(a, n) / norm(cross(a, n));
    circle.ey = cross(n, circle.ex) / norm(cross(n, circle.ex));
end
end
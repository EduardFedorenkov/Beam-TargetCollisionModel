function nuSource = GetNuSource(Nv, vGrid, mg, mp, np, VTp, Vp, diffCrossSection)
nuSource = zeros(Nv*Nv*Nv, Nv*Nv*Nv);
for i = 1:Nv*Nv*Nv
    [ki, li, mi] = ind2sub([Nv, Nv, Nv], i);
    Vi = [vGrid(ki), vGrid(li), vGrid(mi)];
    for j = 1: Nv*Nv*Nv
        [kj, lj, mj] = ind2sub([Nv, Nv, Nv], j);
        Vj = [vGrid(kj), vGrid(lj), vGrid(mj)];

        if (i ~= j)
            nuSource(i, j) = GetSourceFreq(mg, mp, Vi, Vj, np, VTp, Vp, diffCrossSection);
        end
    end
end
end


function nuSource = GetNuSource(Vi_list, normUji_matrix, mg, mp, np, VTp, Vp, diffCrossSection)
Ntotal = size(Vi_list, 1);
nuSource = zeros(Ntotal, Ntotal);
for i = 1:Ntotal
    Vi = Vi_list(i, :);
    nuSourceRow = zeros(1, Ntotal);
    for j = 1: Ntotal
        if i == j
            continue;
        end
        Vj = Vi_list(j, :);
        normUji = normUji_matrix(i, j);

        nuSourceRow(j) = GetSourceFreq(mg, mp, Vi, Vj, normUji, np, VTp, Vp, diffCrossSection);
    end
    nuSource(i, :) = nuSourceRow;
end
end


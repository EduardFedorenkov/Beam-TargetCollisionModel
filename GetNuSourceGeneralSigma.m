function nuSource = GetNuSourceGeneralSigma(Vi_list, normUji_matrix, mg, mp, np, VTp, Vp, Rp)
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

        nuSourceRow(j) = GetSourceFreqGeneralSigma(mg, mp, Vi, Vj, normUji, np, VTp, Vp, Rp);
    end
    nuSource(i, :) = nuSourceRow;
end
end


function nuSink = GetNuSinkMatGeneralSigma(Vi_list, normUji_matrix, mg, mp, np, VTp, Vp, Rp)
Ntotal = size(Vi_list, 1);
nuSink = zeros(Ntotal, Ntotal);
for i = 1:Ntotal
    Vi = Vi_list(i, :);
    nuSinkRow = zeros(1, Ntotal);
    for j = 1: Ntotal
        if i == j
            continue;
        end
        Vj = Vi_list(j, :);
        normUji = normUji_matrix(i, j);

        nuSinkRow(j) = GetSinkFreqGeneralSigma(mg, mp, Vi, Vj, normUji, np, VTp, Vp, Rp);
    end
    nuSink(i, :) = nuSinkRow;
end
end


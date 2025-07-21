function nuSource = GetNuSourceConstSigma(Vi_list, mg, mp, np, VTp, Vp, Rp)
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

        nuSourceRow(j) = GetSourceFreqConstSigma(mg, mp, Vi, Vj, np, VTp, Vp, Rp);
    end
    nuSource(i, :) = nuSourceRow;
end
end


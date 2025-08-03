function nuSink = GetNuSinkMatConstSigma(Vi_list, mg, mp, np, VTp, Vp, Rp)
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

        nuSinkRow(j) = GetSinkFreqConstSigma(mg, mp, Vi, Vj, np, VTp, Vp, Rp);
    end
    nuSink(i, :) = nuSinkRow;
end
end


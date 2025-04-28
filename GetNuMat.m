function [nu1, nu2] = GetNuMat(mg, mp, v_grid, diffCrossSec, np, VTp, Vp, Nv, Nangle)
nu1 = zeros(Nv^3,Nv^3);
nu2 = zeros(Nv^3,Nv^3);

massFactor = (mg + mp) / (2 * mp);
velGridStep = v_grid(2) - v_grid(1);

for kj = 1:Nv
    for ki = 1:Nv
        for lj = 1:Nv
            for li = 1:Nv
                for mj = 1:Nv
                    for mi = 1:Nv
                        idxI = [mi, li, ki];
                        idxJ = [mj, lj, kj];
                        if all(idxI ~= idxJ)
                            Vi = [v_grid(mi), v_grid(li), v_grid(ki)];
                            Vj = [v_grid(mj), v_grid(lj), v_grid(kj)];
                            firstIdx = mi + (li - 1) * Nv + (ki - 1) * Nv * Nv;
                            secondIdx = mj + (lj - 1) * Nv + (kj - 1) * Nv * Nv;
                            nu1(firstIdx, secondIdx) = AngularIntegration1(massFactor, velGridStep, diffCrossSec, Vi, Vj, np, VTp, Vp, Nangle);
                            nu2(firstIdx, secondIdx) = AngularIntegration2(massFactor, velGridStep, diffCrossSec, Vi, Vj, np, VTp, Vp, Nangle);
                        end
                    end
                end
            end
        end
    end
end

end


function I = AngularIntegration1(massFactor, velGridStep, diffCrossSec, Vi, Vj, np, VTp, Vp, Nangle)
    Uji = Vj - Vi;
    UjiNorm = norm(Uji);
    factor = massFactor^2 * UjiNorm * velGridStep^3 * diffCrossSec / ( Nangle / (4 * pi) );
    factor = factor * np / (pi^(3/2) * VTp^3);
    points = RandSampleSphere(Nangle);

    I = 0;
    for i = 1:Nangle
        nz = dot(Uji / UjiNorm, points(i,:));
        if nz > 0
            vel = norm( Vi - Vp + massFactor * ( 2 * Uji - UjiNorm * points(i,:) / nz ) );
            I = I + exp( - vel^2 / VTp^2 ) / nz^3;
        end
    end
    I = I * factor;
end


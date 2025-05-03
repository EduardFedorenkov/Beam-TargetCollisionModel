function fNext = ImplicitTimeScheme(f, nuSource, nuSink, dt)
A = diag(1 + dt * nuSink) - dt * nuSource;
fNext(:) = A \ f(:);
end


function fNext = SemiImplicitTimeScheme(f, nuSource, nuSink, dt)
fNext(:) = (f(:) + dt * nuSource * f(:)) ./ (1 + dt * nuSink);
end


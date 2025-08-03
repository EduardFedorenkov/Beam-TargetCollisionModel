function fNext = ExplicitTimeScheme(f, nuSource, nuSink, dt)
fNext(:) = f(:) + dt * (-nuSink .* f(:) + nuSource * f(:));
end


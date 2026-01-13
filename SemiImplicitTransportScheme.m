function fNext = SemiImplicitTransportScheme(f, nuSource, nuSink, dl, Vp)
fNext(:) = (f(:) + (dl / Vp) * nuSource * f(:)) ./ (1 + (dl / Vp) * nuSink);
end


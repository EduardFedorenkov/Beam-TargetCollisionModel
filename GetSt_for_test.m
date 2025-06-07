function beta = GetSt_for_test(nuSink, nuSource, f)
StSource = nuSource * f(:);
StSink = -nuSink .* f(:);
sourceFull = sum(StSource, 'all');
sinkFull = sum(StSink, 'all');
beta = -sinkFull / sourceFull;
end


function St = GetSt(nuSink, nuSource, f, Nv)
St_tmp = -nuSink .* f(:) + nuSource * f(:);
St = reshape(St_tmp, [Nv, Nv, Nv]);
end


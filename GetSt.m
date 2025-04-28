function St = GetSt(nu1, nu2, f, Nv)
St_tmp = -sum(nu1, 2) .* f(:) + nu2' * f(:);
St = reshape(St_tmp, [Nv, Nv, Nv]);
end


function Ta = TermalizationSolver(ma, mb, Za, Zb, na, nb, Ta0, Tb, dt, N)
    t = 0:dt:(N * dt);
    % Инициализация массива для решения
    n = length(t);
    Ta = zeros(1, n);
    Ta(1) = Ta0;
    
    f = @(T) PlasmaIonsTermolization(ma, mb, Za, Zb, na, nb, T, Tb) * (Tb - T);
    
    % Цикл по временным шагам
    for i = 1:n-1
        Tcurrent = Ta(i);
        
        % Вычисляем коэффициенты Рунге-Кутта
        k1 = f(Tcurrent);
        k2 = f(Tcurrent + dt/2 * k1);
        k3 = f(Tcurrent + dt/2 * k2);
        k4 = f(Tcurrent + dt * k3);
        
        % Обновляем значение T_a
        Ta(i+1) = Tcurrent + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
end


% Физические параметры
c = 2.99792458e10;                                  % Speed of light [cm/s]
eVtoErg = 1.602176634e-12;                          % Convertion coef from [eV] to [Erg]
m = 938.272e6;                                      % Mass of proton [eV]

% Plasma parameters
mp = m;                                             % Ions mass [eV]
np = 1e13;                                          % Ions density [cm^{-3}]
Tp = 3;                                             % Ions temperature [eV]
VTp = sqrt(2 * Tp / mp) * c;                        % Ions termal vel [cm / s]
% VyMean = 0.7 * VTp;
VyMean = 0;

% Wall
Twall = 1;                                          % Wall temperature [eV]
VTwall = sqrt(2 * Twall / mp) * c;                  % Wall termal vel [cm / s]

% Пространственная сетка (2D)
R = 1.0;                                            % Радиус круглой области
Nx = 50; Ny = 50;
x = linspace(-1.2 * R, 1.2 * R, Nx);
y = linspace(-1.2 * R, 1.2 * R, Ny);
[X, Y] = meshgrid(x, y);
dx = x(2) - x(1);

% Сетка по скоростям (3D)
v_max = 5 / sqrt(2) * VTp;                          % Максимальная скорость
Nv = 11;                                            % Число точек по каждой компоненте скорости
vx = linspace(-v_max, v_max, Nv);
vy = linspace(-v_max, v_max, Nv);
vz = linspace(-v_max, v_max, Nv);
dv = vx(2) - vx(1);

dt = 4 * dx / v_max;                                   % Шаг по времени
T = 50 * dt;                                             % Время моделирования

% Начальное условие: Гауссов пучок в центре с разбросом по скоростям
f = zeros(Ny, Nx, Nv, Nv, Nv);
for ivx = 1:Nv
    for ivy = 1:Nv
        for ivz = 1:Nv
            f(:,:,ivx,ivy,ivz) = np / (pi^(3/2) * VTp^3) * exp(-(X.^2 + Y.^2)/0.4^2) * ...
                                 exp(-(vx(ivx)^2 + (vy(ivy) - VyMean)^2 + vz(ivz)^2)/VTp^2);
        end
    end
end

% Предварительное вычисление плотности для проверки нормировки
initial_density = sum(sum(sum(f, 5), 4), 3) * dv^3;
disp(['Initial average density: ', num2str(mean(initial_density(:)))]);

for t = 0:dt:T
    f_new = zeros(size(f));  % Новая функция распределения
    
    for ivx = 1:Nv
        for ivy = 1:Nv
            for ivz = 1:Nv
                vx_curr = vx(ivx);
                vy_curr = vy(ivy);
                
                % Новые координаты (перенос для ВСЕХ частиц)
                X_new = X - vx_curr * dt;
                Y_new = Y - vy_curr * dt;
                
                % Интерполяция новой функции распределения
                f_slice = squeeze(f(:,:,ivx,ivy,ivz));
                f_new_slice = interp2(X, Y, f_slice, X_new, Y_new, 'linear', 0);
                f_new(:,:,ivx,ivy,ivz) = f_new_slice;
                
                % Проверка выхода за границы (только для нарушителей)
                mask_out = (X_new.^2 + Y_new.^2) >= R^2;
                
                if any(mask_out(:))
                    [row, col] = find(mask_out);
                    
                    for p = 1:length(row)
                        x0 = X(row(p), col(p));
                        y0 = Y(row(p), col(p));
                        
                        % Находим точку пересечения с границей
                        a = roots([vx_curr^2 + vy_curr^2, ...
                                  2*(x0*vx_curr + y0*vy_curr), ...
                                  x0^2 + y0^2 - R^2]);
                        a = min(a(a > 0));

                        % Генерируем новую скорость (Максвелл)
                        u1 = rand();
                        u2 = rand();
                        v_new = sqrt(-2 * log(u1)) * cos(2 * pi * u2);
                        v_new = abs(v_new) * VTwall;
                        
                        % Случайное направление (от стенки внутрь)
                        phi = pi / 2 * rand();       % Полярный угол (только полусфера)
                        theta = 2 * pi * rand(); % Азимутальный угол
                        
                        % Новые компоненты скорости
                        vx_new = v_new * sin(phi) * cos(theta);
                        vy_new = v_new * sin(phi) * sin(theta);
                        vz_new = v_new * cos(phi);
                        
                        % Находим ближайшую ячейку сетки по скоростям
                        [~, ivx_new] = min(abs(vx - vx_new));
                        [~, ivy_new] = min(abs(vy - vy_new));
                        [~, ivz_new] = min(abs(vz - vz_new));
                        
                        % Коррекция: убираем вылетевшую частицу и добавляем отраженную
                        f_new(row(p), col(p), ivx, ivy, ivz) = 0;
                        f_new(row(p), col(p), ivx_new, ivy_new, ivz_new) = ...
                            f_new(row(p), col(p), ivx_new, ivy_new, ivz_new) + ...
                            f_slice(row(p), col(p)) * (VTp^3/VTwall^3) * exp(-v_new^2/VTwall^2) / exp(-(vx_curr^2+vy_curr^2+vz(ivz)^2)/VTp^2);
                    end
                end
            end
        end
    end
    
    f = f_new;
    
    % Визуализация
    if mod(t, 5*dt) < dt
        f_density = squeeze(sum(sum(sum(f, 5), 4), 3)) * dv^3;
        figure(1);
        pcolor(x, y, f_density);
        axis equal tight;
        colormap("jet");
        colorbar;
        shading flat;
        shading interp;
        title(['Time = ', num2str(t), ', Average density: ', num2str(mean(f_density(:)))]);
        drawnow;
    end
end

% Проверка сохранения плотности
final_density = sum(sum(sum(f, 5), 4), 3) * dv^3;
disp(['Final average density: ', num2str(mean(final_density(:)))]);


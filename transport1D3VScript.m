% Физические параметры
c = 2.99792458e10;      % Скорость света [см/с]
eV = 1.602176634e-12;   % 1 эВ в эргах
mp = 938.272e6;         % Масса протона [эВ/c^2]

% Параметры задачи
L = 10;                 % Длина области [см]
n0 = 1e13;              % Начальная концентрация [см^-3]
T0 = 2;                 % Температура левой стенки [эВ]
TL = 1;                 % Температура правой стенки [эВ]
VT0 = sqrt(2*T0/mp)*c;  % Термическая скорость (левая стенка)
VTL = sqrt(2*TL/mp)*c;  % Термическая скорость (правая стенка)

% Сетка по пространств3V3у
Nx = 100;
x = linspace(0, L, Nx);
dx = x(2) - x(1);

% Сетка по скоростям (только vx, так как 1D)
Nv = 13;
v_max = 5*max(VT0,VTL);
vx = linspace(-v_max, v_max, Nv);
vy = linspace(-v_max, v_max, Nv);
vz = linspace(-v_max, v_max, Nv);
dv = vx(2) - vx(1);
[Vx, Vy, Vz] = meshgrid(vx, vy, vz);

% Шаг по времени
dt = 0.5*dx/v_max;
T_final = 1000*dt;

% Начальное распределение
f = zeros(Nx, Nv, Nv, Nv);
for ix = 1 : Nx
    f(ix,:,:,:) = n0/(sqrt(pi)*VT0)^3 * exp(-(x(ix) - L/2).^2/3^2) * exp(-(Vx.^2 + Vy.^2 + Vz.^2)/VT0^2);
end

i = 1;

% Основной цикл
for t = 0:dt:T_final
    f_new = zeros(size(f));
    flux_out_left = 0;
    flux_out_right = 0;
    
    % Перенос частиц и учет вылетающих
    for ivx = 1:Nv
        x_new = x - vx(ivx)*dt;
        
        temp_f = reshape(f(:,ivx,:,:), [length(x), Nv*Nv]);
        temp_f_new = interp1(x, temp_f, x_new, 'linear', 0);
        
        % Проверка границ для текущей скорости vx
        if vx(ivx) < 0
            mask_out_left = (x_new <= 0);
            flux_out_left = flux_out_left + sum(temp_f_new(mask_out_left,:).*abs(vx(ivx)), 2)*dv^3;
            temp_f_new(mask_out_left,:) = zeros(length(mask_out_left(mask_out_left~=0)), Nv*Nv);
        elseif vx(ivx) > 0
            mask_out_right = (x_new >= L);
            flux_out_right = flux_out_right + sum(temp_f_new(mask_out_right,:).*vx(ivx), 2)*dv^3;
            temp_f_new(mask_out_right,:) = zeros(length(mask_out_right(mask_out_right~=0)), Nv*Nv);
        end

        f_new(:,ivx,:,:) = reshape(temp_f_new, [length(x_new), 1, Nv, Nv]);
    end
    
    % Добавление отраженных частиц
    if flux_out_left > 0
        % Нормализация для левой стенки
        norm_factor = flux_out_left/( sum(exp(-vx(vx>0).^2/VT0^2)) * ...
            sum(exp(-vy.^2/VT0^2)) * sum(exp(-vz.^2/VT0^2)) * dv^3);
        for ivx = find(vx > 0)
            for ivy = 1:Nv
                for ivz = 1:Nv
                    f_new(1,ivx,ivy,ivz) = f_new(1,ivx,ivy,ivz) + norm_factor * ...
                        exp(-(vx(ivx)^2 + vy(ivy)^2 + vx(ivz)^2)/VT0^2);
                end
            end
        end
    end
    
    if flux_out_right > 0
        % Нормализация для правой стенки
        norm_factor = flux_out_right/( sum(exp(-vx(vx<0).^2/VTL^2)) * ...
            sum(exp(-vy.^2/VTL^2)) * sum(exp(-vz.^2/VTL^2)) * dv^3);
        for iv = find(vx < 0)
            for ivy = 1:Nv
                for ivz = 1:Nv
                    f_new(end,ivx,ivy,ivz) = f_new(end,ivx,ivy,ivz) + norm_factor * ...
                        exp(-(vx(ivx)^2 + vy(ivy)^2 + vx(ivz)^2)/VTL^2);
                end
            end
        end
    end
    
    f = f_new;
    
    % Визуализация
    if mod(t, 100*dt) < dt
        n = sum(sum(sum(f,4),3),2)*dv^3;
        
        figure(i)
        i = i + 1;
        plot(x, n)
        xlabel('x')
        ylabel('Плотность n(x)')
        title(['Время: ' num2str(t)])
    end
end
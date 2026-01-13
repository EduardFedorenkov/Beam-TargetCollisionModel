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

% Сетка по пространств
Nx = 100;
x = linspace(0, L, Nx);
dx = x(2) - x(1);

% Сетка по скоростям (только vx, так как 1D)
Nv = 51;
v_max = 5*max(VT0,VTL);
vx = linspace(-v_max, v_max, Nv);
dv = vx(2) - vx(1);

% Шаг по времени
dt = 0.5*dx/v_max;
T_final = 1000*dt;

% Начальное распределение
f = zeros(Nx, Nv);
for iv = 1:Nv
    f(:,iv) = n0/(sqrt(pi)*VT0) * exp(-(x - L/2).^2/3^2) * exp(-vx(iv)^2/VT0^2);
end

i = 1;

% Основной цикл
for t = 0:dt:T_final
    f_new = zeros(size(f));
    flux_out_left = 0;
    flux_out_right = 0;
    
    % Перенос частиц и учет вылетающих
    for iv = 1:Nv
        x_new = x - vx(iv)*dt;
        f_new(:,iv) = interp1(x, f(:,iv), x_new, 'linear', 0);
        
        % Проверка границ для текущей скорости
        if vx(iv) < 0  % Частицы движутся вправо
            mask_out_left = (x_new <= 0);
            flux_out_left = flux_out_left + sum(f_new(mask_out_left,iv).*abs(vx(iv)))*dv;
            f_new(mask_out_left,iv) = 0;  % Удаляем вылетевшие
        elseif vx(iv) > 0  % Частицы движутся влево
            mask_out_right = (x_new >= L);
            flux_out_right = flux_out_right + sum(f_new(mask_out_right,iv).*vx(iv))*dv;
            f_new(mask_out_right,iv) = 0;  % Удаляем вылетевшие
        end
    end
    
    % Добавление отраженных частиц
    if flux_out_left > 0
        % Нормализация для левой стенки
        norm_factor = flux_out_left/sum(exp(-vx(vx>0).^2/VT0^2)*dv);
        for iv = find(vx > 0)
            f_new(1,iv) = f_new(1,iv) + norm_factor*exp(-vx(iv)^2/VT0^2);
        end
    end
    
    if flux_out_right > 0
        % Нормализация для правой стенки
        norm_factor = flux_out_right/sum(exp(-vx(vx<0).^2/VTL^2)*dv);
        for iv = find(vx < 0)
            f_new(end,iv) = f_new(end,iv) + norm_factor*exp(-vx(iv)^2/VTL^2);
        end
    end
    
    f = f_new;
    
    % Визуализация
    if mod(t, 100*dt) < dt
        n = sum(f,2)*dv;
        J = sum(f*vx',2)*dv;
        
        figure(i)
        i = i + 1;
        subplot(2,1,1)
        plot(x, n)
        xlabel('x')
        ylabel('Плотность n(x)')
        title(['Время: ' num2str(t)])
        
        subplot(2,1,2)
        plot(x, J)
        xlabel('x')
        ylabel('Поток J(x)')
        drawnow
    end
end
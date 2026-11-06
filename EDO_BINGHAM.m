clear; clc; close all; 

%% 
%% ========= EDO ===========
% (1/r) d/dr ( r * tau ) = G ->  d/dr ( r * tau ) = G * r

% ===== Datos =====
mu1 = 0.909;        % ajusta
vartheta   = 680.6; % ajusta
Dp=-64000;          % [Pa] Caida de presion
L=0.3;           
G  = Dp/L;          % ajusta
R   = 7.5e-3;       % 7.5 mm
tau_w=800;
a = R * vartheta / tau_w; % a -> para definir la frontera del comportamiento como fluido de bingham

eps0  = 1e-9;                  % evitar r=0 exacto (singularidad)
rmesh = linspace(a, R, 80); % Malla de trabajo de la funcion bvp4c (EDO)

%% ===== PREPARACION DE LA EDO =====
% Para poder trabajar con la EDO debemos pasar la ecuacion de 2 orden a 1 orden

% y = [v; dvdr] -> "y" es la derivada de "v"
% y(1)= y1(r)= v
% y(2)= y2(r)= v'=dv/dr
edo_fun = @(r,y) [ y(2);
                     ((G/mu1)*r + (vartheta/mu1) - y(2))./r ]; %Recibe dos entradas llamadas "r"(radio en nuestro caso) e "y" (vector de estado y(1) y(2))

%% ===== CONDICIONES DE CONTORNO DE LA EDO =====
% BC: dv/dr(0)=0, v(R)=0
% ya: vector y evaluado en r=a
% yb: vector y evaluado en r=R
% Podrian denominarse como CONDICIONES DE CONTORNO[v'(a)=0; v(R)=0 No-Slip]
fun_contorno = @(ya,yb) [ ya(2);
                          yb(1) ];

%% ===== SOLUCION DE LA EDO =====
solinit = bvpinit(rmesh, [0;0]); % rmesh: primeros puntos tentativos de resolucion. [0,0] conjetura de partida.
sol = bvp4c(edo_fun, fun_contorno, solinit);


%% ===== Construir la solución completa (tapón + región cortante)=====
% Evaluación en una malla completa (0..R)
rplot = linspace(eps0, R, 300);

% Velocidad en r=a. OJO: sol está definido en [a, R]
va = deval(sol, a);
va = va(1);   % v(a)

vplot = zeros(size(rplot));

% Tapón: r <= a -> v constante = v(a)
plug = (rplot <= a);
vplot(plug) = va;

% Región cortante: r > a -> evaluar la solución numérica
no_plug = (rplot > a);
y_no_plug = deval(sol, rplot(no_plug));   % devuelve [v; dv/dr]
vplot(no_plug) = y_no_plug(1,:);

% Simetría respecto a r = 0
rplot_sym = [-fliplr(rplot), rplot];
vplot_sym = [ fliplr(vplot), vplot];

plot(rplot_sym*1e3, vplot_sym, 'LineWidth', 2), grid on
xlabel('r [mm]')
ylabel('v_z(r)')
title('Perfil de velocidad con núcleo no deformado (Bingham)')

%% ========= Tabla de velocidades en radios específicos (mm) =========
% Radios de interés (mm)
r_mm = [-7.5;-7.4;-7.3;-7.2;-7.1;-7;-6.8;-6.6;-6.4;-6.2;-6;-5.5;-5;-4.5;-4;-3.5;-3;-2.5;-2;-1.5;-1;-0.5; 0;0.5;1;1.5;2;2.5;3;3.5;4;4.5;5;5.5;6;6.2;6.4;6.6;6.8;7;7.1;7.2;7.3;7.4;7.5];

% Convertimos a metros
r_m = r_mm * 1e-3;

% Vector de velocidades
v_z = zeros(size(r_m));

% Interpolador del perfil calculado (solo r >= 0)
v_interp = @(rr) interp1(rplot, vplot, rr, 'pchip');

for i = 1:numel(r_m)
    rr = abs(r_m(i));   % simetría: v(-r) = v(r)

    if rr <= a
        % Zona tapón
        v_z(i) = va;    % velocidad constante del núcleo
    else
        % Zona cedida
        v_z(i) = v_interp(rr);
    end
end

% Crear tabla
T = table(v_z, ...
          'VariableNames', {'v_z'});

disp('Tabla de velocidades v_z en los radios solicitados:')
disp(T)
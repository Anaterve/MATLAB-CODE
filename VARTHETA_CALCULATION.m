%% Bingham laminar en tubo: calcular vartheta (= tau0) a partir de m_dot y dP/dz
% Entradas (define aquí tus datos):
m_dot = 0.0212;      % [kg/s]  caudal masico
Dp=-64000;       % [Pa] Caida de presion 
L = 0.300;      % [m]  Long de la tuberia
G  = Dp/L;     % [Pa/m]  gradiente de presion (normalmente negativo)
R     = 0.0075;      % [m]     radio del tubo
mu    = 0.909;      % [Pa·s]  viscosidad plastica (mu1)
rho   = 1805;      % [kg/m^3] densidad

% Opcional para el perfil:
Nr = 200;          % puntos radiales para v(r)

%% 1) Caudal volumetrico medido
Q = m_dot / rho;   % [m^3/s]

%% 2) Esfuerzo cortante en pared a partir de G
% Para flujo totalmente desarrollado en tubo:
% tau(r) = (r/2) * G  => tauw = (R/2)*G
tauw = abs(G) * R / 2;   % [Pa] (magnitud)

%% 3) Comprobaciones de factibilidad [Mediante Hagen-Poisseuille] (si tau0=0, recupera Newtoniano)
Q_newton = pi * R^4 * abs(G) / (8 * mu);   % [m^3/s] limite superior para tau0=0

if Q <= 0
    error('Q_meas <= 0: revisa m_dot y rho.');
end

if Q > Q_newton
    error(['No hay solucion Bingham (laminar) con estos datos: ', 'Q_meas > Q_newton. Revisa mu, R o G.']);
end

%% 4) Resolver x = vartheta/tauw usando la ecuacion de caudal de Bingham
% QB = (pi*R^3*tauw/(12*mu)) * [ 3 - 4x + x^4 ],  con x=tau0/tauw
K = (pi * R^3 * tauw) / (12 * mu);

f = @(x) K * (3 - 4*x + x.^4) - Q;

% x debe estar en [0, 1). En x->1, QB->0.
x_lo = 0.0;
x_hi = 1.0 - 1e-12;

% Asegurar cambio de signo (por robustez)
f_lo = f(x_lo);     % = Q_newton - Q_meas >= 0
f_hi = f(x_hi);     % ~ -Q < 0

if ~(f_lo >= 0 && f_hi < 0)
    error('No se puede acotar la raiz en [0,1). Revisa datos.');
end

% Resolver con fzero usando un intervalo (bracketing)
x = fzero(f, [x_lo, x_hi]);

%% 5) Calcular vartheta y radio del nucleo (plug)
vartheta = x * tauw;     % tauw [Pa]
a = R * x;               % [m] radio del plug

%% 6) Perfil de velocidad v_z(r) segun tu formulacion (opcional)
r = linspace(0, R, Nr).';
vz = zeros(size(r));

% Region de corte: R > r > a
idxShear = (r > a) & (r <= R);
rr = r(idxShear);

vz(idxShear) = (tauw * R / (2*mu)) .* (1 - (rr./R).^2) ...
             - (vartheta * R / mu) .* (1 - (rr./R));

% Region plug: r <= a (constante)
vz_plug = (tauw * R / (2*mu)) * (1 - (vartheta/tauw))^2;
vz(r <= a) = vz_plug;

% Velocidad media a partir del caudal medido
Vmean = Q / (pi * R^2);

%% 7) Mostrar resultados
fprintf('--- Resultados Bingham laminar ---\n');
fprintf('Q_meas      = %.6e m^3/s\n', Q);
fprintf('tauw        = %.6f Pa\n', tauw);
fprintf('vartheta(tau0) = %.6f Pa\n', vartheta);
fprintf('x=tau0/tauw  = %.6f [-]\n', x);
fprintf('a (radio plug) = %.6e m\n', a);
fprintf('Vmean       = %.6e m/s\n', Vmean);

%% 8) Plot del perfil (opcional)
figure; plot(r, vz); grid on
xlabel('r [m]'); ylabel('v_z(r) [m/s]');
title('Perfil de velocidad Bingham laminar en tubo');
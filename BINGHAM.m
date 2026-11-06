clear; clc; close all; 

%% 1. DATOS (pon aquí tus valores)
% Gval = dp/dz < 0 (caída de presión a lo largo de +z)
Gval   = -64000/0.3;    % G = dp/dz
vartheta  = 680.6;        % Esfuerzo de cedencia (vartheta)
muval  = 0.909;     % viscosidad plástica
Rval   = 0.0075;    % radio del conducto

%% 2. VARIABLES SIMBÓLICAS
syms r G th mu R C1 C2 a
syms tau(r) v(r)

%% ========= EDO 1 =========================
% (1/r) d/dr ( r * tau ) = G ->  d/dr ( r * tau ) = G * r

% Integramos: r*tau = ∫ G*r dr + C1
int1   = int(G*r, r);          % = G*r^2/2
r_tau  = int1 + C1;            % r*tau = G*r^2/2 + C1

% tau(r) antes de imponer regularidad
tau_rz = simplify(r_tau / r); % = G*r/2 + C1/r

% Regularidad en el eje (r->0): eliminamos el término C1/r
tau_rz_sol = subs(tau_rz, C1, 0);  % tau(r) = G*r/2  (G<0 -> tau<0)

disp('EDO 1 resuelta: tau_rz(r) ='); pretty(tau_rz_sol)

%% Radio de cedencia a: tau(a) = -vartheta  ->  tau(a) = r/2*G
eq_a   = subs(tau_rz_sol, r, a) == -th;   % G*a/2 = -th
aval = solve(eq_a, a);               % a = -2*th/G  (positivo porque G<0)

disp('a simbólico ='); pretty(aval)

%% ========= EDO 2: Bingham =============
% En a <= r <= R, sabemos que dv/dr < 0, tau<0.
% Modelo Bingham (signo coherente con dv/dr<0):
%   tau_rz = -vartheta + mu * dv/dr
%   mu * dv/dr = tau_rz + vartheta
%   dv/dr      = (tau_sol + th)/mu -> Tensor de BINGHAM

dvdr_expr = (tau_rz_sol + th)/mu;

% Integramos: v(r) = ∫ dv/dr dr + C2
v_int = int(dvdr_expr, r) + C2;

% Condición de contorno: no deslizamiento en la pared v(R) = 0
eq_C2   = subs(v_int, r, R) == 0;
C2_expr = solve(eq_C2, C2);

% Sustituimos C2 -> solución v(r) zona cedencia
v_shear_expr = simplify( subs(v_int, C2, C2_expr) );

disp('EDO 2 resuelta: v(r) en zona cedencia =');
pretty(v_shear_expr)

%% Zona tapón (0 <= r <= a): dv/dr = 0 -> v = constante = v(a)
v_plug_expr = simplify( subs(v_shear_expr, r, aval) );
disp('v_plug (cte) ='); pretty(v_plug_expr)

%% ========= Sustitución de valores numéricos =======================

Gsym  = sym(Gval);
thsym = sym(vartheta);
musym = sym(muval);
Rsym  = sym(Rval);

a_num = double( subs(aval, [G, th], [Gsym, thsym]) );

v_shear_num(r) = subs(v_shear_expr, [G, th, mu, R], [Gsym, thsym, musym, Rsym]);

v_plug_num = double( subs(v_plug_expr, [G, th, mu, R], [Gsym, thsym, musym, Rsym]) );

v_shear_fun = matlabFunction(v_shear_num, 'Vars', r);

fprintf('Radio de cedencia a = %.6g m\n', a_num);
fprintf('Velocidad del tapón v_plug = %.6g\n', v_plug_num);

%% ========= Perfil en 0 <= r <= R (dos tramos) =====================

% Tramo tapón: 0 <= r <= a
N1 = 100;
r_plug = linspace(0, a_num, N1);
v_plug = v_plug_num * ones(size(r_plug));

% Tramo cedido: a <= r <= R
N2 = 200;
r_shear = linspace(a_num, Rval, N2);
v_shear_vals = v_shear_fun(r_shear);

% Juntamos r >= 0
r_pos = [r_plug, r_shear];
v_pos = [v_plug, v_shear_vals];

%% ========= Extensión simétrica [-R, R] y gráfica ==================

% Simetría radial: v(-r) = v(r)
r_all = [-fliplr(r_pos(2:end)), r_pos];   % no repetimos r=0
v_all = [ fliplr(v_pos(2:end)), v_pos];

figure;
plot(v_all, r_all, 'LineWidth', 2); hold on;
yline(a_num,  '--k', ' a');
yline(-a_num, '--k', '-a');
grid on;
xlabel('v_z(r)');
ylabel('r');
title('Fluido Bingham: EDO esfuerzos + EDO velocidad (tapón y zona cedida)');
ylim([-Rval Rval]);

% Save the plot as a PNG file
% saveas(gcf, 'Bingham_Fluid_Profile.png');

%% ========= Tabla de velocidades en radios específicos (mm) =========

% Radios de la foto (mm)
r_mm = [-7.4;-6;-4.5;-3;-1.5;0;1.5;3;4.5;6;7.4];

% Convertimos a metros 
r_m = r_mm * 1e-3;

% Perfil de velocidad con simetría: v(-r)=v(r)
v_tab = zeros(size(r_m));
for i = 1:numel(r_m)
    rr = abs(r_m(i));               % simetría radial
    if rr <= a_num                  % zona tapón
        v_tab(i) = v_plug_num;
    else                            % zona cedida
        v_tab(i) = v_shear_fun(rr);
    end
end

% Creamos y mostramos tabla
%T = table(r_mm, r_m, v_tab, 'VariableNames', {'r_mm','r_m','v_z'});
T = table(v_tab, 'VariableNames', {'v_z'})
disp('Tabla de velocidades v_z en los radios solicitados:');
disp(T);  
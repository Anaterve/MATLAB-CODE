%% 
%% ==============================================================
%  RESOLUCION NUMERICA IMPLICITA 1D MEDIANTE DIFERENCIAS FINITAS
%
%  Ecuacion:
%  d/dr [ r ( mu dvz/dr - theta(r,t) ) ] = r G(r,t) - r rho dvz/dt
%
%  Se emplea:
%    - Diferencias finitas en la variable radial r
%    - Metodo de Euler implicito en el tiempo
%    - Condicion de simetria en r = 0
%    - Condicion de Dirichlet en r = R. Fijar un valor de la función en r=R
%    (en nuestro caso, la condicion de no deslizamiento)
%
% ==============================================================

clear; clc; close all;

%% ==============================================================
%  1. PARAMETROS DEL PROBLEMA
% ==============================================================
rho = 1805;     % densidad
mu  = 0.909;    % viscosidad
theta_y = theta_dato;

R = 0.0075;      % radio maximo
T = 1.0;        % tiempo final

N  = 100;       % numero de subintervalos espaciales en direccion radial
Nt = 5000;       % numero de pasos temporales

dr = R / N;     % paso espacial radial. Espacio entre dos nodos consecutivos.
dt = T / Nt;    % paso temporal. Tiempo entre dos instantes consecutivos.

r = linspace(0, R, N+1)';      % Vector columna con todos los nodos radiales
t = linspace(0, T, Nt+1);      % Vector fila para todos los nodos temporales

%% ==============================================================
%  2. RESERVA DE MEMORIA
% ==============================================================
% V(i,n) aproxima v_z(r_i, t^n)
V = zeros(N+1, Nt+1);    % Matriz donde se guarda la velocidad axial     

%% ==============================================================
%  3. CONDICION INICIAL
% ==============================================================
for i = 1:N+1
    V(i,1) = condicion_inicial(r(i));    % Aqui se rellena la primera columna de V con la condicion final definida al final del codigo
end

%% ==============================================================
%  4. BUCLE TEMPORAL
% ==============================================================
for n = 1:Nt %En cada iteracion temporal se calcula la solucion en el instante de t(n+1) usando la solucion conocida en t(n)
    
    tn1 = t(n+1);   % tiempo t^(n+1). 
    
     
    A = zeros(N+1, N+1);    % Matriz del sistema lineal
    b = zeros(N+1, 1);    %Termino independiente
    % Contiene:
    %      - informacion del tiempo anterior V(:,n)
    %      - terminos fuente G(r,t)
    %      - terminos conocidos de theta(r,t)
    %      - condiciones de contorno
   
    % ----------------------------------------------------------
    % 4.1 CONDICION DE SIMETRIA EN r = 0
    %
    % dv/dr = 0  =>  V(1,n+1) = V(2,n+1)
    % Esta ecuacion se coloca en la primera fila de A.
    % ----------------------------------------------------------
    A(1,1) = 1.0;
    A(1,2) = -1.0;
    b(1)   = 0.0;
    
    % ----------------------------------------------------------
    % 4.2 NODOS INTERIORES
    %
    % Los nodos interiores son los que no estan ni en el centro r = 0 ni en la pared r = R.
    % Para j = 2,...,N corresponden a los nodos fisicos interiores.
    % En cada nodo se discretiza la ecuacion diferencial original.
    % ----------------------------------------------------------
    for j = 2:N 
        
        ri  = r(j);        % Posicion radial del nodo actual
        rip = ri + dr/2;   % representa r_(i+1/2). Es el punto medio entre el nodo actual y el nodo de la derecha.
        rim = ri - dr/2;   % representa r_(i-1/2). Es el punto medio entre el nodo actual y el nodo de la izquierda.
        
        %% ------------------------------------------------------
        %  Evaluacion de theta en las caras del volumen de control
        % ------------------------------------------------------
        % La ecuacion tiene el termino: d/dr [ r * theta(r,t) ]
        % Por eso se necesita theta en los puntos intermedios:
        %      theta_{i+1/2}
        %      theta_{i-1/2}
        %
        % Como theta se conoce en los nodos, se aproxima en las caras
        % usando un promedio aritmetico.
        % Valores medios de theta en t^(n+1)
        m_reg = 100;   % parametro de regularizacion

        dvdr_ip = (V(j+1,n) - V(j,n)) / dr;
        dvdr_im = (V(j,n)   - V(j-1,n)) / dr;
        
        mu_app_ip = mu + theta_y * (1 - exp(-m_reg*abs(dvdr_ip))) / ...
                         (abs(dvdr_ip) + 1e-12);
        
        mu_app_im = mu + theta_y * (1 - exp(-m_reg*abs(dvdr_im))) / ...
                 (abs(dvdr_im) + 1e-12);
        
        % Coeficientes de la matriz A
        % En cada nodo interior, la ecuacion discreta tiene la forma:
        %      aW * V(j-1,n+1) + aP * V(j,n+1) + aE * V(j+1,n+1) = rhs
        % donde:
        %      aW = coeficiente del nodo de la izquierda
        %      aP = coeficiente del nodo central
        %      aE = coeficiente del nodo de la derecha
        %
        % Estos coeficientes vienen de discretizar el termino viscoso:
        %
        %      d/dr [ r * mu * dvz/dr ]
        %
        % y el termino transitorio:
        %
        %      r * rho * dvz/dt

        aW = - mu_app_im * rim / dr^2;   % Coeficiente asociado al nodo izquierdo. 
        % Es negativo porque aparece como parte de la segunda derivada radial discretizada.
        
        aP = ri * rho / dt ...
   + (mu_app_ip*rip + mu_app_im*rim) / dr^2;
        % Coeficiente asociado al nodo central.
        % Tiene dos contribuciones:
        % 1) ri * rho / dt - Procede del termino transitorio implicito.
        % 2) mu * (rip + rim) / dr^2 - Procede del termino viscoso.
        
        aE = - mu_app_ip * rip / dr^2;   % Coeficiente asociado al nodo derecho.
        % Tambien es negativo por la discretizacion del termino viscoso.
        
        
        % Termino independiente. Contiene todo lo conocido en el instante actual.
                rhs = (ri * rho / dt) * V(j,n) ...
            + ri * G_dato(ri, tn1);
        % Primer termino:
        %     (ri * rho / dt) * V(j,n)
        % Es la parte conocida del termino temporal.
        % Euler implicito:
        %      dvz/dt ≈ [V^(n+1) - V^n] / dt
        % La parte con V^(n+1) va a la matriz A.
        % La parte con V^n queda en rhs.

        % Segundo termino:
        %      ri * G_dato(ri, tn1)
        % Representa la fuerza fuente G evaluada en el nodo actual y en el tiempo nuevo.

        % Tercer termino:
        %      (rip * theta_ip - rim * theta_im) / dr
        % Aproxima:
        %      d/dr [ r * theta(r,t) ]
        % usando diferencias finitas en las caras del volumen radial.

               
        % Ensamblaje de la fila j
        A(j,j-1) = aW;  % Coloca el coeficiente del nodo izquierdo.
        A(j,j)   = aP;  % Coloca el coeficiente del nodo central.
        A(j,j+1) = aE;  % Coloca el coeficiente del nodo derecho.
        b(j)     = rhs; % Coloca el termino independiente correspondiente a la fila j.
    end
    
    % ----------------------------------------------------------
    % 4.3 CONDICION DE CONTORNO EN r = R
    %
    % v(R,t) = h(t)
    % ----------------------------------------------------------
    % En la pared de la tuberia se impone condicion de no deslizamiento.
    % Es decir, la velocidad del fluido en contacto con la pared es igual a la velocidad de la pared.
    % Si la pared esta quieta:
    %      v_z(R,t) = 0
    % Esta condicion se impone sustituyendo la ultima ecuacion del sistema por:
    %      V(N+1,n+1) = condicion_borde(tn1)
    
    A(N+1, N+1) = 1.0;
    b(N+1)      = condicion_borde(tn1);
    
    % ----------------------------------------------------------
    % 4.4 RESOLUCION DEL SISTEMA LINEAL
    % ----------------------------------------------------------
   
    % Se resuelve:
    %      A * V(:,n+1) = b
    % El operador "\" de MATLAB resuelve sistemas lineales.
    % El resultado es el perfil radial completo en el tiempo nuevo.

    V(:,n+1) = A \ b;
end

%% ==============================================================
%  5. REPRESENTACION GRAFICA
% ==============================================================

% Perfil final
figure;
plot(r, V(:,end), 'LineWidth', 2);
xlabel('r');
ylabel('v_z(r,T)');
title('Perfil radial final');
grid on;

% Evolucion espacio-tiempo
figure;
surf(t, r, V, 'EdgeColor', 'none');
xlabel('t');
ylabel('r');
zlabel('v_z(r,t)');
title('Solucion numerica mediante esquema implicito');
colorbar;
view(45,30);



%% ==============================================================
%  6. FUNCIONES DEL PROBLEMA
%
%  Sustituir por los datos concretos del enunciado
% ==============================================================
    % 6.1 VALOR DE VARTHETA
    % 6.1.1 Constante
    function val = theta_dato(r,t)
    % Funcion theta(r,t) conocida
    %val = 0.2 * r .* exp(-t);
    val = 0;
    end
    


    % 6.2 VALOR DE G
    % 6.2.1 Entrada en rampa. G solo depende de t no de r.
    function val = G_dato(r,t)
    % Funcion G(r,t) con entrada en rampa semibrusca
    tr = 0.1;      % tiempo de subida de la rampa. Durante los primeros 0.1 s, G aumenta desde 0 hasta A
    A  = 19500.0;   % valor final de G (Pa/m)

    if t < 0
        g_t = 0;
    elseif t <= tr
        g_t = A * (t/tr);
    else
        g_t = A;
    end

    val =g_t;
    end
    
    % 6.2.2 Cambio suave / presion creciente
    %function val = G_dato(r,t)
    %R   = 0.075;   % radio maximo
    %G0  = 50.0;     % amplitud final
    %tau = 0.2;     % tiempo caracteristico de transicion
    %val = G0 * (1 - exp(-t/tau));
    %end

    % 6.2.3 Caida de presion
    %function val = G_dato(r,t)
    %R   = 0.075;    % radio maximo
    %G0  = 50.0;     % amplitud final
    %tc = 0.1;       % instante en que empieza a desaparecer
    %tau = 0.03;     % suavidad de la trnsicion
    %val = G0 * 0.5 * (1 - tanh((t - tc)/tau));
    %end
    
    % 6.2.4 Cambio en el valor de P por variacion en las cond. del flujo
    % function val = G_dato(r,t)
    % G0 = 19500.0;     % P del Poiseuille inicial
    % G1 = 225500;      % Nuevo valor de Presion
    % 
    % if t <= 0
    %     val = G0;
    % else
    %     val = G1;
    % end
    % end


    % 6.3 VALOR INICIAL DE v_z
    % 6.3.1 Flujo parado (arranque)
    function val = condicion_inicial(r)
    % Condicion inicial v_z(r,0)
    val = 0.0 * r;
    end
    
    % 6.3.2 Poisseuille
    % function val = condicion_inicial(r)
    % % Perfil de Poiseuille con v_max = 0.3 m/s
    % R = 0.0075;        % radio maximo
    % v_max = 0.3;
    % val = v_max * (1 - (r./R).^2);
    % end



    % 6.3 CONDICION DE BORDE
    % 6.3.1 No deslizamiento
    function val = condicion_borde(t)
    % Condicion de contorno en r = R
    val = 0.0;
    end
    
    %% ==============================================================
    %  7. ESTABILIZACION DEL SISTEMA
    %  Se consideró que el sistema alcanzaba el régimen estacionario cuando 
    % la diferencia relativa entre dos perfiles consecutivos era inferior a 10−4
    % en norma infinita.
    % ==============================================================

%% ==============================================================
%  7.1. ANALISIS DE CONVERGENCIA AL ESTADO ESTACIONARIO
% ==============================================================

% Perfil estacionario aproximado
V_final = V(:,end);

% Error relativo respecto al perfil final
error_rel = zeros(1,Nt+1);

for n = 1:Nt+1

    error_rel(n) = ...
        max(abs(V(:,n) - V_final)) ...
        / max(abs(V_final));

end

% Grafica de convergencia
figure;
semilogy(t,error_rel,'LineWidth',2);
xlabel('Tiempo (s)');
ylabel('Error relativo');
title('Convergencia al estado estacionario');
grid on;

% Tiempo de estabilizacion
tol = 1e-3;      % tolerancia del 0.1 %

idx = find(error_rel < tol,1,'first');

if ~isempty(idx)

    t_estacionario = t(idx);

    fprintf('\n');
    fprintf('====================================\n');
    fprintf('TIEMPO DE ESTABILIZACION (perfil)\n');
    fprintf('t_est = %.6f s\n',t_estacionario);
    fprintf('====================================\n');

else

    fprintf('\n');
    fprintf('La simulacion no alcanza la tolerancia.\n');

end


%% ==============================================================
%  7.2. ANALISIS DE CONVERGENCIA ENTRE PASOS TEMPORALES
% ==============================================================

error_dt = zeros(1,Nt);

for n = 2:Nt+1

    error_dt(n-1) = ...
        max(abs(V(:,n)-V(:,n-1))) ...
        / max(abs(V(:,n)));

end

% Grafica
figure;
semilogy(t(2:end),error_dt,'LineWidth',2);
xlabel('Tiempo (s)');
ylabel('Cambio relativo');
title('Variacion entre pasos temporales');
grid on;

% Tiempo de estabilizacion
tol_dt = 1e-6;

idx_dt = find(error_dt < tol_dt,1,'first');

if ~isempty(idx_dt)

    t_estacionario_dt = t(idx_dt+1);

    fprintf('\n');
    fprintf('====================================\n');
    fprintf('TIEMPO DE ESTABILIZACION (incremental)\n');
    fprintf('t_est = %.6f s\n',t_estacionario_dt);
    fprintf('====================================\n');

else

    fprintf('\n');
    fprintf('La simulacion no alcanza la tolerancia incremental.\n');

end
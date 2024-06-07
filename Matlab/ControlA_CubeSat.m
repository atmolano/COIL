%% Control Actitud CubeSat

clear all
clc

%% Constantes, coeficientes y matrices

% Constantes
r_orbita = 6378.14e3 + 500e3; % Radio órbita CubeSat (m)
G_tierra = 398600.442e9; % Constante gravitacional de la Tierra (m^3/s^2)
w_orbita = sqrt(G_tierra / r_orbita^3); % Frecuencia angular de la órbita (rad/s)

% Momentos de inercia del CubeSat (kg*m^2) 

m_inercia_X = 1.391e-6;
m_inercia_Y = 1.467e-6;
m_inercia_Z = 1.565e-6;

% Coeficientes de la matriz A
coef1 = 3 * w_orbita^2 * (m_inercia_Z - m_inercia_Y) / m_inercia_X;
coef2 = w_orbita * (m_inercia_Y - m_inercia_Z) / m_inercia_X;
coef3 = 2 * w_orbita^2 * (m_inercia_Z - m_inercia_X) / m_inercia_Y;
coef4 = w_orbita * (m_inercia_X - m_inercia_Y) / m_inercia_Z;

% Matriz A para el sistema dinámico del CubeSat
matrizA = [0            0            -w_orbita  1            0   0;
           0            0            0          0            1   0;
           w_orbita     0            0          0            0   1;
           coef1 0            0          0            0   coef2;
           0            coef3 0          0            0   0;
           0            0            0          coef4 0   0];

% Matriz de entradas B
matrizB = [0                0                0;
           0                0                0;
           0                0                0;
           1/m_inercia_X    0                0;
           0                1/m_inercia_Y    0;
           0                0                1/m_inercia_Z];

% Campo magnético de la Tierra
campoMag_Tierra = [10; 15; 20];

matrizBFinal = zeros(size(matrizB));

% Multiplicación por componente
for i = 1:3
    matrizBFinal(:, i) = matrizB(:, i) * campoMag_Tierra(i);
end

% Calculo polos 
p = [-1.5 -1.6 -1.7 -1.8 -1.9 -2.0];

% Valores de la diagonal
diag_vals = [0.8, 0.9, 1.1];

%matriz B nueva
matriz_B_new = zeros(size(matrizBFinal));

% Multiplicación elemento por elemento
for i = 1:numel(diag_vals)
    matriz_B_new(:, i) = matrizBFinal(:, i) * diag_vals(i);
end


% Ganancias del controlador usando los polos
GanananciasControlador = place(matrizA, matriz_B_new, p);

%% Condiciones iniciales y solución de ecuaciones diferenciales

% Condiciones iniciales  CubeSat
ang_R0 = 0.3; % radianes
ang_P0 = -0.1; % radianes
ang_Y0 = 0.6; % radianes
w_R0 = 0.02; % radianes/s
w_P0 = -0.01; % radianes/s
w_Y0 = 0.03; % radianes/s

% Solución ecuaciones diferenciales
[tiempo, estados] = ode45(@(tiempo, estado) Eq_Mov_CubeSat(tiempo, estado, GanananciasControlador), [0 100], [ang_R0 ang_P0 ang_Y0 w_R0 w_P0 w_Y0]);

% Extracción de ángulos y tasas angulares
h_ang   = estados(:, 1:3);
h_w = estados(:, 4:6);

%% Gráficas 

% Grafica ángulos
figure
subplot(2, 1, 1)
plot(tiempo, h_ang)
title('Ángulos de Actitud')
xlabel('Tiempo (s)')
ylabel('Ángulo (rad)')
legend('Ángulo R (\phi)', 'Ángulo P (\theta)', 'Ángulo Y (\psi)')
grid on

% Grafica tasas angulares
subplot(2, 1, 2)
plot(tiempo, h_w)
title('Velocidades Angulares')
xlabel('Tiempo (s)')
ylabel('Velocidad (rad/s)')
legend('Velocidad R (\omega_x)', 'Velocidad P (\omega_y)', 'Velocidad Y (\omega_z)')
grid on

%% Funciones

% Función para las ecuaciones de movimiento del CubeSat
function derivadas = Eq_Mov_CubeSat(tiempo, estado, GananciasControlador)
    % Momentos de inercia del CubeSat (kg*m^2)
    m_inercia_X = 1.391e-6;
    m_inercia_Y = 1.467e-6;
    m_inercia_Z = 1.565e-6;

    % Descomposición de estado
    [ang_R, ang_P, ang_Y, w_R, w_P, w_Y] = deal(estado(1), estado(2), estado(3), estado(4), estado(5), estado(6));

    % Controlador
    in_control = -GananciasControlador * estado;

    % Ecuaciones de movimiento
    dev_Ang_R = w_R;
    dev_Ang_P = w_P;
    dev_Ang_Y = w_Y;
    dev_Vel_R = (m_inercia_Y - m_inercia_Z) * w_P * w_Y / m_inercia_X + in_control(1) / m_inercia_X;
    dev_Vel_P = (m_inercia_Z - m_inercia_X) * w_R * w_Y / m_inercia_Y + in_control(2) / m_inercia_Y;
    dev_Vel_Y = (m_inercia_X - m_inercia_Y) * w_R * w_P / m_inercia_Z + in_control(3) / m_inercia_Z;

    % Derivadas
    derivadas = [dev_Ang_R; dev_Ang_P; dev_Ang_Y; dev_Vel_R; dev_Vel_P; dev_Vel_Y];
end



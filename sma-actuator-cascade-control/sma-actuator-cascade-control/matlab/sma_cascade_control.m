%% SMA ACTUATOR - CASCADE CONTROL DESIGN
% Position control of a load lifted by a Shape Memory Alloy (SMA) wire
% through a pulley. Inner loop: wire temperature (voltage -> T).
% Outer loop: load position (T -> z).
%
% Run this script first: it computes the equilibrium, the linearised plants
% and the two controllers, and defines the workspace variables used by
% sma_nonlinear_closed_loop.slx (A1..D1 inner controller, A2..D2 outer
% controller, equilibrium values, physical parameters).

clear
close all
clc

%% Figure settings
font_size   = 15;
line_width  = 2;
plot_pos    = [500, 300, 600, 400];

s = tf('s');

%% Simulink solver settings
MaxStep = 1e-5;   % at least 10x smaller than the fastest time constant
RelTol  = 1e-5;
t_f     = 10;     % simulation stop time [s]

%% Electrical parameters
Re = 4.8;         % wire electrical resistance [Ohm]
Le = 55e-3;       % inductance [H]

%% Thermal parameters
T0   = 273.15;    % 0 degC in K
c    = 6.2;       % steepness of the stiffness transition
l    = 15e-2;     % wire length [m]
d    = 0.2e-3;    % wire diameter [m]
Cth  = 25.1e-3;   % thermal capacity [J/K]
Rth  = 16.7;      % thermal resistance towards the Peltier cell [K/W]
Tp   = 15 + T0;   % Peltier cell temperature [K]
h    = 6.55;      % convection coefficient [W/(m^2 K)]
Tamb = 24 + T0;   % ambient temperature [K]

%% Mechanical parameters
r1   = 0.45e-2;   % inner pulley radius [m]
r2   = 4.5e-2;    % outer pulley radius [m]
J    = 2e-4;      % pulley inertia [kg m^2]
m    = 15e-3;     % load mass [kg]
beta = 1.5;       % viscous friction [N s/m]

%% SMA stiffness K(T) = Kmax*(1 - 1/(1 + exp(c*(T - Tavg)/Tdiff)))
Kmax  = 3.92e3;   % [N/m]
Tavg  = 70;       % [degC]
Tdiff = 20;       % [degC]
x3_eq = 4e-2;     % desired equilibrium position z* [m]

%% Derived parameters
m_eq = J/r2 + m*r2;
g    = 9.81;
Ramb = 1/(h*pi*d*l);         % equivalent convection thermal resistance
S    = pi*d*l;               % lateral wire surface
rho  = r2/r1;                % pulley ratio
deltal_max = 4/100*l;        % maximum wire contraction

%% Equilibrium point
K_star  = (m*g*r2)/(r1*(deltal_max - x3_eq*r1/r2));
x2_eq   = Tavg + (Tdiff/c)*log(K_star/(Kmax - K_star));          % T* [degC]
x2_eq_k = x2_eq + T0;                                            % T* [K]
x1_eq   = sqrt((1/Re)*(h*S*(x2_eq_k - Tamb) + (1/Rth)*(x2_eq_k - Tp)));  % I*
x4_eq   = 0;
u_eq    = Re*x1_eq;                                              % V*

fprintf('Equilibrium: I* = %.4f A, T* = %.2f degC, z* = %.3f m, V* = %.4f V\n', ...
        x1_eq, x2_eq, x3_eq, u_eq);

%% Full linearised model (state x = [I, T, z, dz/dt], input V, output z)
expo   = exp(c*(x2_eq_k - (Tavg + T0))/(Tdiff + T0));
K_diff = (c*Kmax/(Tdiff + T0))*(expo/((1 + expo)^2));
K_x2eq = Kmax*(1 - 1/(1 + expo));

A = [-Re/Le            0                      0                                 0;
     2*Re*x1_eq/Cth   -1/Cth*(h*S + 1/Rth)    0                                 0;
     0                 0                      0                                 1;
     0   ((r1*(deltal_max - x3_eq*r1/r2))/(r2*m_eq))*K_diff   -((r1^2)/((r2^2)*m_eq))*K_x2eq   -beta/m_eq];
B = [1/Le; 0; 0; 0];
C = [0 0 1 0];
D = 0;
G = tf(ss(A, B, C, D));

%% Specifications
% - zero steady-state error for step reference (0.5 cm) and step dTp (5 degC)
% - phase margin > 65 deg
% - attenuation > 800x of sinusoidal dTp (5 degC) below 5 rad/s
% - attenuation > 10x of sinusoidal measurement noise (0.05 cm) above 20 rad/s
% - overshoot < 5 %, 1 % settling time < 1.8 s
Ta = 1.7;         % design settling time [s]

%% Inner-loop plant Gt(s): voltage -> temperature (linearised at I*)
dAt = [-Re/Le             0;
        2*Re*x1_eq/Cth   -1/(Rth*Cth) - 1/(Ramb*Cth)];
dBt = [1/Le   0;
       0      1/(Rth*Cth)];
dCt = [0 1];
dDt = 0;
Gt = tf(ss(dAt, dBt, dCt, dDt));   % Gt(1): V -> T,  Gt(2): dTp -> T
Gt(1)

%% Outer-loop plant Gm(s): temperature -> position (linearised at T*, z*)
n = exp(c/Tdiff*(x2_eq - Tavg));
dAm = [0                                          1;
       -Kmax*(1 - 1/(1 + n))*r1/(m_eq*rho)        -beta*r2/m_eq];
dBm = [0;
       Kmax*r1/m_eq*(deltal_max - x3_eq/rho)*c/Tdiff*n/(1 + n)^2];
dCm = [1 0];
dDm = 0;
Gm = tf(ss(dAm, dBm, dCm, dDm))

%% Inner-loop controller: PID  R1 = mu (1 + tauz s)^2 / (s (1 + taup s)(1 + 1.2e-5 s))
% The design places the double zero to obtain the target margin Mf1 at wc1.
% The two high-frequency poles added afterwards cost about 6 deg of phase,
% so the target is set above the 65 deg specification.
% (Exam version: Mf1 = 65 -> actual margin ~59 deg.)
wc1 = 290;        % desired crossover [rad/s]
Mf1 = 72;         % design phase margin [deg]

[M1, P1, W1] = bode(Gt(1)/s, wc1);
ModGt = 20*log10(M1);
ArgGt = P1;
if ArgGt > 0, ArgGt = ArgGt - 360; end   % keep the unwrapped phase
Pdi   = -180 + Mf1 - ArgGt;           % phase lead required from the two zeros
rho_star1 = tand(Pdi/2);
tauz  = rho_star1/wc1;
A1dB  = 20*log10(sqrt(1 + rho_star1^2));
u1    = 10^(-(2*A1dB + ModGt)/20);
taup  = 1/(10*wc1);
R1 = u1*(1 + tauz*s)^2/s/(1 + taup*s)/(1 + 1.2e-5*s);
L1 = R1*Gt(1);

F1 = minreal(L1/(1 + L1));            % complementary sensitivity
S1 = minreal(1/(1 + L1));             % sensitivity
Q1 = minreal(R1/(1 + L1));            % control sensitivity
R1

%% Outer-loop controller: PI  R2 = mu (1 + tau s) / (s (1 + 0.024 s)(1 + 1e-4 s))
% Designed on Gm with the inner loop treated as ideal (F1 ~ 1 at these
% frequencies). Target raised for the same reason as above
% (exam version: Mf = 65 -> actual margin ~59 deg, overshoot ~9 %).
Mf = 75;
wc = 460/(Ta*65);                     % crossover from the settling-time spec

[M2, P2, W2] = bode(Gm/s, wc);
ModGm = 20*log10(M2);
ArgGm = P2;
if ArgGm > 0, ArgGm = ArgGm - 360; end
Pdo   = -180 + Mf - ArgGm;
rho_star = tand(Pdo);
tau   = rho_star/wc;
Ao    = 20*log10(sqrt(1 + rho_star^2));
u     = 10^(-(Ao + ModGm)/20);
R2 = u*(1 + tau*s)/(s*(1 + 0.024*s)*(1 + 0.0001*s));
L2 = R2*Gm;

F2 = minreal(R2*F1*Gm/(1 + R2*F1*Gm));
S2 = minreal(1/(1 + R2*F1*Gm));
Q2 = minreal(R2/(1 + R2*F1*Gm));
R2

%% Margins and step-response metrics
[~, Pm1, ~, Wc1] = margin(L1);
[~, Pm2, ~, Wc2] = margin(R2*F1*Gm);
info2 = stepinfo(F2, 'SettlingTimeThreshold', 0.01);
fprintf('Inner loop: wc = %.0f rad/s, PM = %.1f deg\n', Wc1, Pm1);
fprintf('Outer loop: wc = %.2f rad/s, PM = %.1f deg, overshoot = %.1f %%, Ts(1%%) = %.2f s\n', ...
        Wc2, Pm2, info2.Overshoot, info2.SettlingTime);

%% State-space realisations used by the Simulink model
[A1, B1, C1, D1] = ssdata(ss(R1));    % inner controller
[A2, B2, C2, D2] = ssdata(ss(R2));    % outer controller

%% Plots
% Control effort, inner loop
figure(1); set(1, 'position', plot_pos)
bodemag(Gt(1), R1, Q1); grid on
set(findobj(1, 'Type', 'Line'), 'LineWidth', line_width)
title('Control sensitivity (inner loop)', 'FontSize', font_size)
legend('G_T', 'R_1', 'Q_1', 'FontSize', font_size)

% Control effort, outer loop
figure(2); set(2, 'position', plot_pos)
bodemag(Gm, R2, Q2); grid on
set(findobj(2, 'Type', 'Line'), 'LineWidth', line_width)
title('Control sensitivity (outer loop)', 'FontSize', font_size)
legend('G_M', 'R_2', 'Q_2', 'FontSize', font_size)

% Bode of L1 with the dTp attenuation constraint (> 58 dB below 5 rad/s)
w = logspace(-1, 4, 1000);
[mag, phase] = bode(L1, w);
figure(3)
subplot(2, 1, 1)
semilogx(w, 20*log10(squeeze(mag))); hold on
patch([0.1 0.1 5 5], [60 -60 -60 60], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
title('Bode diagram - inner loop L_1'); ylabel('Magnitude [dB]'); grid on
subplot(2, 1, 2)
semilogx(w, squeeze(phase)); hold on
xlabel('Frequency [rad/s]'); ylabel('Phase [deg]'); grid on

% Bode of L2 with settling-time and noise constraints
w = logspace(-1, 2, 1000);
[mag, phase] = bode(L2, w);
figure(4)
subplot(2, 1, 1)
semilogx(w, 20*log10(squeeze(mag))); hold on
patch([0.1 0.1 4 4], [0 -60 -60 0], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
patch([20 20 100 100], [60 -20 -20 60], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
title('Bode diagram - outer loop L_2'); ylabel('Magnitude [dB]'); grid on
subplot(2, 1, 2)
semilogx(w, squeeze(phase)); hold on
patch([4 4 20 20], [-115 -180 -180 -115], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Frequency [rad/s]'); ylabel('Phase [deg]'); grid on

% Sensitivity functions
figure(5); set(5, 'position', plot_pos)
bodemag(F1, S1, Q1); grid on
set(findobj(5, 'Type', 'Line'), 'LineWidth', line_width)
title('Sensitivity functions (inner loop)', 'FontSize', font_size)
legend('F_1', 'S_1', 'Q_1', 'FontSize', font_size)

figure(6); set(6, 'position', plot_pos)
bodemag(F2, S2, Q2); grid on
set(findobj(6, 'Type', 'Line'), 'LineWidth', line_width)
title('Sensitivity functions (outer loop)', 'FontSize', font_size)
legend('F_2', 'S_2', 'Q_2', 'FontSize', font_size)

% Step responses
figure(7); set(7, 'position', plot_pos)
[y, t] = step(F1);
plot(t, y, 'LineWidth', line_width); grid on; xlim([0, t(end)])
xlabel('Time [s]', 'FontSize', font_size); ylabel('Temperature [K]', 'FontSize', font_size)
title('Closed-loop step response (inner loop)', 'FontSize', font_size)

figure(8); set(8, 'position', plot_pos)
[y, t] = step(F2);
plot(t, y, 'LineWidth', line_width); grid on; xlim([0, t(end)])
xlabel('Time [s]', 'FontSize', font_size); ylabel('Position [m]', 'FontSize', font_size)
title('Closed-loop step response (outer loop)', 'FontSize', font_size)

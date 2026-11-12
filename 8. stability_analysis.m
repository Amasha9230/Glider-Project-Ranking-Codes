%% LONGITUDINAL STATIC STABILITY, DOWNWASH & CG DETERMINATION PIPELINE
% Free-Flight Glider Project | Department of Mechanical Engineering, UoM
clear; clc; close all;

%% =========================================================================
% 1. LOCKED DESIGN INPUTS (From Wing & Tail Finalization Stages)
% =========================================================================
% Main Wing Parameters (Drela AG14, AR = 6.80)
b_w      = 0.400;       % Wingspan (m)
S_w      = 0.02353;     % Reference Wing Area (m^2)
AR_w     = 6.80;        % Wing Aspect Ratio
c_root   = 0.0713;      % Wing Root Chord (m)
c_tip    = 0.0464;      % Wing Tip Chord (m)
c_bar    = 0.0597;      % Mean Aerodynamic Chord (MAC) (m)
lambda_w = 0.65;        % Wing Taper Ratio
x_ac_w   = 0.01493;     % Wing AC location behind MAC leading edge (m) (0.25 * c_bar)

% Horizontal Tail Parameters (Drela HT12, Sized via flow5 Ranking)
S_t      = 0.00370;     % Horizontal Tail Area (m^2)
AR_t     = 3.50;        % Tail Aspect Ratio
l_t      = 0.2100;      % Tail Moment Arm (from wing AC to tail AC) (m)
c_t_bar  = 0.0325;      % Tail Mean Chord (m)
x_ac_t   = 0.25 * c_t_bar; % Tail AC location behind tail leading edge (m)

% Flight Baseline & Target Static Margin
V_cruise = 6.0;         % Cruise Velocity (m/s)
rho      = 1.225;       % Air Density (kg/m^3)
SM_des   = 0.10;        % Target Static Margin (10% MAC)

%% =========================================================================
% 2. FORMULATION & CALCULATION ENGINE
% =========================================================================
fprintf('================================================================================\n');
fprintf('       LONGITUDINAL STATIC STABILITY & BALANCE CALCULATION ENGINE                \n');
fprintf('================================================================================\n\n');

% -------------------------------------------------------------------------
% STEP A: 3D Finite-Wing Lift Curve Slopes (Helmbold Formulation)
% -------------------------------------------------------------------------
% 2D theoretical sectional lift slopes for low-Re profiles (rad^-1)
% Drela AG14 ~ 5.85 rad^-1 | Drela HT12 ~ 5.908 rad^-1
a0_w = 5.850;   % Wing 2D sectional lift slope (rad^-1)
a0_t = 5.908;   % Tail 2D sectional lift slope (rad^-1)

% 3D finite wing lift curve slope: CL_alpha = a0 / (1 + (57.3 * a0) / (pi * e * AR))
% Using Helmbold equation: a = a0 / [sqrt(1 + (a0 / (pi * AR))^2) + (a0 / (pi * AR))]
CL_alpha_w = a0_w / (sqrt(1 + (a0_w / (pi * AR_w))^2) + (a0_w / (pi * AR_w)));
CL_alpha_t = a0_t / (sqrt(1 + (a0_t / (pi * AR_t))^2) + (a0_t / (pi * AR_t)));

CL_alpha_w_deg = CL_alpha_w * (pi / 180); % per degree
CL_alpha_t_deg = CL_alpha_t * (pi / 180); % per degree

fprintf('[STEP 1] 3D FINITE-WING LIFT-CURVE SLOPES:\n');
fprintf('  Formula : CL_alpha = a0 / [sqrt(1 + (a0 / (pi * AR))^2) + (a0 / (pi * AR))]\n');
fprintf('  Inputs  : Wing a0 = %.3f rad^-1, AR_w = %.2f | Tail a0 = %.3f rad^-1, AR_t = %.2f\n', ...
    a0_w, AR_w, a0_t, AR_t);
fprintf('  Outputs : Wing CL_alpha = %.4f rad^-1 (%.4f deg^-1)\n', CL_alpha_w, CL_alpha_w_deg);
fprintf('          : Tail CL_alpha = %.4f rad^-1 (%.4f deg^-1)\n\n', CL_alpha_t, CL_alpha_t_deg);

% -------------------------------------------------------------------------
% STEP B: Wing Downwash Gradient at Tail Station (d_epsilon / d_alpha)
% -------------------------------------------------------------------------
% Classical low-speed downwash formulation (Nelson / Roskam):
% d_eps / d_alpha = (2 * CL_alpha_w) / (pi * AR_w)
deps_dalpha = (2 * CL_alpha_w) / (pi * AR_w);
downwash_factor = 1 - deps_dalpha;

fprintf('[STEP 2] DOWNWASH GRADIENT AT EMPENNAGE STATION:\n');
fprintf('  Formula : d_epsilon / d_alpha = (2 * CL_alpha_w) / (pi * AR_w)\n');
fprintf('  Formula : Downwash Factor     = (1 - d_epsilon / d_alpha)\n');
fprintf('  Inputs  : CL_alpha_w = %.4f rad^-1, AR_w = %.2f\n', CL_alpha_w, AR_w);
fprintf('  Outputs : d_epsilon / d_alpha = %.4f\n', deps_dalpha);
fprintf('          : Downwash Factor     = %.4f (Dynamic pressure flow factor)\n\n', downwash_factor);

% -------------------------------------------------------------------------
% STEP C: Tail Volume Coefficient (V_h)
% -------------------------------------------------------------------------
V_h = (S_t * l_t) / (S_w * c_bar);

fprintf('[STEP 3] HORIZONTAL TAIL VOLUME COEFFICIENT:\n');
fprintf('  Formula : V_h = (S_t * l_t) / (S_w * c_bar)\n');
fprintf('  Inputs  : S_t = %.5f m^2, l_t = %.4f m, S_w = %.5f m^2, c_bar = %.4f m\n', ...
    S_t, l_t, S_w, c_bar);
fprintf('  Output  : V_h = %.4f (Target stable band: 0.40 - 0.60)\n\n', V_h);

% -------------------------------------------------------------------------
% STEP D: Aircraft Neutral Point Location (x_NP)
% -------------------------------------------------------------------------
% Total aircraft neutral point behind MAC leading edge:
% x_NP = x_ac_w + [(CL_alpha_t / CL_alpha_w) * (1 - deps_dalpha) * (S_t / S_w) * l_t]
tail_contribution = (CL_alpha_t / CL_alpha_w) * downwash_factor * (S_t / S_w) * l_t;
x_NP = x_ac_w + tail_contribution;
x_NP_fraction = (x_NP / c_bar) * 100; % % of MAC

fprintf('[STEP 4] AIRCRAFT NEUTRAL POINT (x_NP) CALCULATION:\n');
fprintf('  Formula : x_NP = x_ac_w + (CL_alpha_t / CL_alpha_w) * (1 - d_eps/d_alpha) * (S_t / S_w) * l_t\n');
fprintf('  Inputs  : x_ac_w = %.4f m (%.2f mm), S_t/S_w = %.4f, l_t = %.4f m\n', ...
    x_ac_w, x_ac_w * 1000, S_t / S_w, l_t);
fprintf('  Outputs : Tail Arm Shift Delta_x = +%.4f m (+%.2f mm)\n', ...
    tail_contribution, tail_contribution * 1000);
fprintf('          : Neutral Point (x_NP)   = %.4f m (%.2f mm behind MAC LE)\n', ...
    x_NP, x_NP * 1000);
fprintf('          : Non-Dimensional x_NP   = %.2f%% MAC\n\n', x_NP_fraction);

% -------------------------------------------------------------------------
% STEP E: Center of Gravity Location for 10% Static Margin (x_CG)
% -------------------------------------------------------------------------
% SM = (x_NP - x_CG) / c_bar = 0.10 => x_CG = x_NP - (SM_des * c_bar)
x_CG = x_NP - (SM_des * c_bar);
x_CG_fraction = (x_CG / c_bar) * 100; % % of MAC

fprintf('[STEP 5] REQUIRED CENTER OF GRAVITY (x_CG) BALANCING:\n');
fprintf('  Formula : x_CG = x_NP - (Static_Margin * c_bar)\n');
fprintf('  Inputs  : Static Margin = %.2f (10%%), c_bar = %.4f m (%.2f mm)\n', ...
    SM_des, c_bar, c_bar * 1000);
fprintf('  Outputs : Required CG Position x_CG = %.4f m (%.2f mm behind MAC LE)\n', ...
    x_CG, x_CG * 1000);
fprintf('          : Non-Dimensional x_CG     = %.2f%% MAC\n\n', x_CG_fraction);

% -------------------------------------------------------------------------
% STEP F: Restoring Pitch Moment Derivative (Cm_alpha)
% -------------------------------------------------------------------------
% Total aircraft pitching moment slope:
% Cm_alpha = CL_alpha_w * ((x_CG - x_ac_w) / c_bar) - CL_alpha_t * downwash_factor * V_h
Cm_alpha_rad = CL_alpha_w * ((x_CG - x_ac_w) / c_bar) - CL_alpha_t * downwash_factor * V_h;
Cm_alpha_deg = Cm_alpha_rad * (pi / 180);

fprintf('[STEP 6] PITCH RESTORING MOMENT DERIVATIVE (Cm_alpha):\n');
fprintf('  Formula : Cm_alpha = CL_alpha_w * [(x_CG - x_ac_w) / c_bar] - CL_alpha_t * (1 - d_eps/d_alpha) * V_h\n');
fprintf('  Outputs : Cm_alpha = %.4f rad^-1 (%.4f deg^-1)\n', Cm_alpha_rad, Cm_alpha_deg);
if Cm_alpha_rad < 0
    fprintf('  Status  : STRONGLY STATICALLY STABLE (Inherent Nose-Down Restoring)\n');
else
    fprintf('  Status  : UNSTABLE (Divergent)\n');
end
fprintf('================================================================================\n\n');

%% =========================================================================
% 3. AUTOMATED VISUALIZATION: STABILITY DIAGRAM & BALANCE MAP
% =========================================================================
fig = figure('Name', 'Longitudinal Stability & Balance Map', 'Color', 'w', 'Position', [100, 100, 900, 480]);

% Subplot 1: Cm vs Alpha across CG envelope
subplot(1, 2, 1); hold on; grid on; box on;
alpha_range = -4:0.5:8; % deg
alpha_rad   = alpha_range * (pi / 180);

% Trim pitch moment with wing Cm0 = -0.034, decalage = +2 deg
Cm0_total = -0.034 + CL_alpha_t * downwash_factor * V_h * (2.0 * pi / 180);
Cm_curve  = Cm0_total + Cm_alpha_rad * alpha_rad;

plot(alpha_range, Cm_curve, 'b-', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('10%% SM (Cm_\\alpha = %.4f deg^{-1})', Cm_alpha_deg));
yline(0, 'k--', 'LineWidth', 1.0, 'DisplayName', 'Trim (Cm = 0)');
xline(0, 'k:', 'LineWidth', 0.8);

xlabel('Angle of Attack \alpha (deg)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Pitching Moment Coefficient C_m', 'FontSize', 10, 'FontWeight', 'bold');
title('Total Glider Pitch Stability Curve', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 8);

% Subplot 2: Physical Chord Balance Map
subplot(1, 2, 2); hold on; grid on; box on;
barh([1, 2, 3], [x_ac_w, x_CG, x_NP] * 1000, 0.45, 'FaceColor', [0.2 0.5 0.8]);
set(gca, 'YTick', [1, 2, 3], 'YTickLabel', ...
    {'Wing AC (x_{ac,w})', 'Target CG (x_{CG})', 'Neutral Point (x_{NP})'}, ...
    'FontSize', 10, 'FontWeight', 'bold');

xlabel('Position behind MAC Leading Edge (mm)', 'FontSize', 10, 'FontWeight', 'bold');
title(sprintf('Longitudinal Balance Allocation (MAC = %.1f mm)', c_bar * 1000), ...
    'FontSize', 11, 'FontWeight', 'bold');
xlim([0, c_bar * 1000]);

% Annotate Static Margin gap
text(x_CG * 1000 + 0.5, 2.5, sprintf('<- 10%% SM Gap (%.1f mm) ->', (x_NP - x_CG) * 1000), ...
    'Color', 'r', 'FontWeight', 'bold', 'FontSize', 9);

saveas(fig, 'Figure_Static_Stability_Balance.png');
fprintf('[OK] Balance & Stability plot exported as "Figure_Static_Stability_Balance.png"\n\n');
clear; clc; close all;

% Parameters

a = 0.3;                                 % Allee threshold
c_star = sqrt(2) * (0.5 - a);            % critical speed
epsilon = 1e-3;                          % small parameter
c = c_star - epsilon;

F = @(U) U .* (1 - U) .* (U - a);

%%
% 1) One-term TW solution U0(z), V0(z)=dU/dz

[z_U0, U0_sol, V0_sol] = solve_U0(F, c_star);

% Ensure U decreases from 1 to 0
if U0_sol(1) < U0_sol(end)
    z_U0   = flipud(z_U0(:));
    U0_sol = flipud(U0_sol(:));
    V0_sol = flipud(V0_sol(:));
else
    z_U0   = z_U0(:);
    U0_sol = U0_sol(:);
    V0_sol = V0_sol(:);
end

%% 
% 2) Numerical phase-plane solution V(U)

[U_num, V_num] = solve_combined_V(F, c_star, epsilon);

U_num = U_num(:);
V_num = V_num(:);

% Ensure U decreases from 1 to 0
if U_num(1) < U_num(end)
    U_num = flipud(U_num);
    V_num = flipud(V_num);
end

% Keep original V for phase-plane plotting
V_num_phase = V_num;

% Only use a clipped version for z integration to avoid division by zero
V_num_z = V_num;
small_idx = abs(V_num_z) < 1e-10;
V_num_z(small_idx) = -1e-10;   % keep negative sign

% z(U) = int dU / V
z_num = cumtrapz(U_num, 1 ./ V_num_z);

% Align numerical curve with one-term curve at U = 0.5
z_U0_half  = interp1(U0_sol, z_U0, 0.5, 'linear', 'extrap');
z_num_half = interp1(U_num,  z_num, 0.5, 'linear', 'extrap');
z_num = z_num + (z_U0_half - z_num_half);

%%
% 3) Prepare phase-plane ordering

[U_num_pp, idx_num_pp] = sort(U_num, 'ascend');
V_num_pp = V_num_phase(idx_num_pp);

[U0_pp, idx_U0_pp] = sort(U0_sol, 'ascend');
V0_pp = V0_sol(idx_U0_pp);


%%
% 4) Figure (b): U versus z
% Only plot Numerical solution and One-term solution

shift_num = -7.76073;
shift_one = -7.76073;

fig1 = figure('Color','w','Position',[100 100 620 460]);
ax1 = axes(fig1);
hold(ax1,'on');

% Numerical solution: blue solid
h_num = plot(ax1, z_num + shift_num, U_num, ...
    'b-', 'LineWidth', 1.6);

% One-term solution: red dashed
h_one = plot(ax1, z_U0 + shift_one, U0_sol, ...
    'r--', 'LineWidth', 1.6);

xlabel(ax1, '$z$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax1, '$U$', 'Interpreter', 'latex', 'FontSize', 16);

legend(ax1, [h_num, h_one], ...
    {'Numerical solution', 'One-term solution'}, ...
    'Interpreter', 'latex', 'Location', 'northeast', ...
    'FontSize', 9, 'Box', 'on');

xlim(ax1, [-16.5, 1]);
ylim(ax1, [-0.02, 1.02]);
grid(ax1, 'on');
box(ax1, 'on');
set(ax1, 'FontSize', 13, 'LineWidth', 1.0, ...
    'TickLabelInterpreter','latex');

% inset
ax1_inset = axes('Position', [0.24 0.28 0.22 0.22]);
hold(ax1_inset, 'on');

plot(ax1_inset, z_num + shift_num, U_num, ...
    'b-', 'LineWidth', 1.2);

plot(ax1_inset, z_U0 + shift_one, U0_sol, ...
    'r--', 'LineWidth', 1.2);

xlim(ax1_inset, [-16.2, 0.2]);
ylim(ax1_inset, [0, 0.03]);
grid(ax1_inset, 'on');
box(ax1_inset, 'on');
set(ax1_inset, 'FontSize', 8, 'LineWidth', 0.8, ...
    'TickLabelInterpreter','latex');

ax1_inset.XTick = [-16 -14 -12 -10 0];
ax1_inset.YTick = [0 0.01 0.02 0.03];

%%
% 5) Figure (c): phase plane V versus U
% Only plot Numerical solution and One-term solution

fig2 = figure('Color','w','Position',[180 120 620 460]);
ax2 = axes(fig2);
hold(ax2,'on');

% Numerical solution: blue solid
h2_num = plot(ax2, U_num_pp, V_num_pp, ...
    'b-', 'LineWidth', 1.6);

% One-term solution: red dashed
h2_one = plot(ax2, U0_pp, V0_pp, ...
    'r--', 'LineWidth', 1.6);

xlabel(ax2, '$U$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax2, '$V$', 'Interpreter', 'latex', 'FontSize', 16);

legend(ax2, [h2_num, h2_one], ...
    {'Numerical solution', 'One-term solution'}, ...
    'Interpreter', 'latex', 'Location', 'northeast', ...
    'FontSize', 9, 'Box', 'on');

xlim(ax2, [0, 1]);
ylim(ax2, [-0.35, 0.01]);
grid(ax2, 'on');
box(ax2, 'on');
set(ax2, 'FontSize', 13, 'LineWidth', 1.0, ...
    'TickLabelInterpreter','latex');

% inset for small-U region
ax2_inset = axes('Position', [0.18 0.26 0.22 0.22]);
hold(ax2_inset, 'on');

plot(ax2_inset, U_num_pp, V_num_pp, ...
    'b-', 'LineWidth', 1.2);

plot(ax2_inset, U0_pp, V0_pp, ...
    'r--', 'LineWidth', 1.2);

xlim(ax2_inset, [0, 0.02]);
ylim(ax2_inset, [-0.03, 0.01]);
grid(ax2_inset, 'on');
box(ax2_inset, 'on');
set(ax2_inset, 'FontSize', 8, 'LineWidth', 0.8, ...
    'TickLabelInterpreter','latex');

ax2_inset.XTick = [0.005 0.01 0.015 0.02];
ax2_inset.YTick = [-0.03 -0.02 -0.01 0 0.01];

% optional vector export
% exportgraphics(fig1,'fig9b_two_curves.pdf','ContentType','vector');
% exportgraphics(fig2,'fig9c_two_curves.pdf','ContentType','vector');


function [z_sol, U0_sol, V0_sol] = solve_U0(F, c)
    % Solve one-term travelling wave ODE in z:
    % U' = V,  V' = -cV - F(U)

    odefun = @(z, Y) [Y(2); -c*Y(2) - F(Y(1))];

    z_span = linspace(-12, 50, 2500);
    Y0 = [1; -1e-4];
    opts = odeset('RelTol',1e-9,'AbsTol',1e-12,'MaxStep',5e-3);

    [z_sol, Y_sol] = ode15s(odefun, z_span, Y0, opts);

    U0_sol = Y_sol(:,1);
    V0_sol = Y_sol(:,2);

    % Trim unphysical part if overshoot occurs
    keep = isfinite(U0_sol) & isfinite(V0_sol) & ...
           (U0_sol >= -0.05) & (U0_sol <= 1.05);

    z_sol  = z_sol(keep);
    U0_sol = U0_sol(keep);
    V0_sol = V0_sol(keep);
end

function [U_sol, V_combined] = solve_combined_V(F, c, epsilon)
    % Solve phase-plane equation:
    %
    % dV/dU = (-cV - F(U))/V
    %
    % and first-order correction:
    %
    % V = V0 + epsilon V1

    U_values = linspace(1, 1e-6, 2500).';

    dV0dU = @(U, V0) (-c*V0 - F(U)) ./ V0;

    V0_start = -epsilon;   % negative branch
    opts = odeset('RelTol',1e-9,'AbsTol',1e-12,'MaxStep',5e-3);

    [U_sol_V0, V0_sol] = ode15s(dV0dU, U_values, V0_start, opts);

    U_sol_V0 = U_sol_V0(:);
    V0_sol   = V0_sol(:);

    % Remove invalid points
    keep0 = isfinite(U_sol_V0) & isfinite(V0_sol);
    U_sol_V0 = U_sol_V0(keep0);
    V0_sol   = V0_sol(keep0);

    V0_interp = @(Uq) interp1(U_sol_V0, V0_sol, Uq, ...
        'linear', 'extrap');

    dV1dU = @(U, V1) 1 - ...
        ((-c*V0_interp(U) - F(U))./V0_interp(U).^2 + ...
        c./V0_interp(U)) .* V1;

    V1_start = 0;
    [U_sol_V1, V1_sol] = ode15s(dV1dU, U_values, V1_start, opts);

    U_sol_V1 = U_sol_V1(:);
    V1_sol   = V1_sol(:);

    % Interpolate V1 onto V0 grid
    V1_interp = interp1(U_sol_V1, V1_sol, U_sol_V0, ...
        'linear', 'extrap');

    V_combined = V0_sol + epsilon * V1_interp;
    U_sol = U_sol_V0;

    keep = isfinite(U_sol) & isfinite(V_combined);
    U_sol = U_sol(keep);
    V_combined = V_combined(keep);
end
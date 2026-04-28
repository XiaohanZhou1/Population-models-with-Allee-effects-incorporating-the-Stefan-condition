clear; clc; close all;

% Parameters
a = 0.3;                                 % Allee threshold
c_star = sqrt(2) * (0.5 - a);            % critical speed
epsilon = 1e-3;                          % small parameter
c = c_star - epsilon;

F = @(U) U .* (1 - U) .* (U - a);

% Matching point
U_match = 0.05;

%% =========================
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

V_num_phase = V_num;

% z(U) = int dU / V
V_num_z = V_num;
V_num_z(abs(V_num_z) < 1e-10) = -1e-10;
z_num = cumtrapz(U_num, 1 ./ V_num_z);

% Align numerical curve with one-term curve at U = 0.5
z_U0_half  = interp_monotone(U0_sol, z_U0, 0.5);
z_num_half = interp_monotone(U_num,  z_num, 0.5);
z_num = z_num + (z_U0_half - z_num_half);

%%
% 3) Explicit inner small-U asymptotic branch
%
% psi_0(phi) ~ -gamma * phi + D * phi^(-r)
%
% U = epsilon^p * phi
% V = epsilon^q * psi
%
% where
% p = q = r = 1/(2a+1),
% gamma = 1/sqrt(2),
% D = -Gamma(2a+1) Gamma(3-2a) / 6.

p = 1 / (2*a + 1);
q = p;
r = p;

gamma_inner = 1 / sqrt(2);
D = -gamma(2*a + 1) * gamma(3 - 2*a) / 6;

epsp = epsilon^p;
epsq = epsilon^q;

psi0 = @(phi) -gamma_inner .* phi + D .* phi.^(-r);

% Convert matching point to phi variable
phi_match = U_match / epsp;

% Avoid phi = 0 because phi^(-r) is singular
U_floor = 1e-8;
phi_floor = U_floor / epsp;

if phi_floor >= phi_match
    phi_floor = 1e-4 * phi_match;
    U_floor = epsp * phi_floor;
end

N_inner = 3000;

% Inner branch in descending U order: U_match -> U_floor
phi_inner = linspace(phi_match, phi_floor, N_inner).';
U_small = epsp * phi_inner;
V_small = epsq * psi0(phi_inner);

% Remove invalid values
valid_idx = isfinite(U_small) & isfinite(V_small) & (U_small > 0) & (V_small < -1e-12);

U_small = U_small(valid_idx);
V_small = V_small(valid_idx);

if numel(U_small) < 2
    error('Too few valid inner points. Check U_match, U_floor, or the explicit inner approximation.');
end

% Ensure descending U
if U_small(1) < U_small(end)
    U_small = flipud(U_small);
    V_small = flipud(V_small);
end

fprintf('\nInner approximation parameters:\n');
fprintf('p                  = %.8f\n', p);
fprintf('q                  = %.8f\n', q);
fprintf('r                  = %.8f\n', r);
fprintf('gamma              = %.8f\n', gamma_inner);
fprintf('D                  = %.8e\n', D);
fprintf('epsilon^p          = %.8e\n', epsp);
fprintf('epsilon^q          = %.8e\n', epsq);
fprintf('U_match            = %.8e\n', U_match);
fprintf('phi_match          = %.8e\n', phi_match);
fprintf('U_floor            = %.8e\n', U_floor);
fprintf('phi_floor          = %.8e\n', phi_floor);
fprintf('Number inner pts   = %d\n', numel(U_small));
fprintf('U_small range      = [%.8e, %.8e]\n', min(U_small), max(U_small));
fprintf('V_small range      = [%.8e, %.8e]\n', min(V_small), max(V_small));

%%
% 4) Build two-term solution
%
% U >= U_match : one-term solution
% U <  U_match : explicit inner branch
% z coordinate for inner branch

V_small_z = V_small;
V_small_z(abs(V_small_z) < 1e-12) = -1e-12;

z_small = cumtrapz(U_small, 1 ./ V_small_z);

% Values at U_match
V_one_match   = interp_monotone(U0_sol,  V0_sol,  U_match);
z_one_match   = interp_monotone(U0_sol,  z_U0,    U_match);

V_inner_match = interp_monotone(U_small, V_small, U_match);
z_inner_match = interp_monotone(U_small, z_small, U_match);

% Align inner z with one-term z at U_match
z_small = z_small + (z_one_match - z_inner_match);

% outer part: one-term branch
idx_outer = U0_sol > U_match;

U_two_outer = [U0_sol(idx_outer); U_match];
V_two_outer = [V0_sol(idx_outer); V_one_match];
z_two_outer = [z_U0(idx_outer);   z_one_match];

% inner part: explicit inner asymptotic branch
idx_inner = U_small < U_match;

U_two_inner = [U_match; U_small(idx_inner)];
V_two_inner = [V_inner_match; V_small(idx_inner)];
z_two_inner = [z_one_match; z_small(idx_inner)];

%%
% 5) Prepare phase-plane ordering

[U_num_pp, idx_num_pp] = sort(U_num, 'ascend');
V_num_pp = V_num_phase(idx_num_pp);

[U0_pp, idx_U0_pp] = sort(U0_sol, 'ascend');
V0_pp = V0_sol(idx_U0_pp);

[U_two_outer_pp, idx_outer_pp] = sort(U_two_outer, 'ascend');
V_two_outer_pp = V_two_outer(idx_outer_pp);

[U_two_inner_pp, idx_inner_pp] = sort(U_two_inner, 'ascend');
V_two_inner_pp = V_two_inner(idx_inner_pp);

%%
% 6) Diagnostics

fprintf('\nDiagnostics:\n');
fprintf('epsilon            = %.6e\n', epsilon);
fprintf('U_match            = %.6e\n', U_match);
fprintf('Number inner pts   = %d\n', numel(U_two_inner));
fprintf('U_inner range      = [%.6e, %.6e]\n', min(U_two_inner), max(U_two_inner));
fprintf('V_inner range      = [%.6e, %.6e]\n', min(V_two_inner), max(V_two_inner));
fprintf('V_one(U_match)     = %.8e\n', V_one_match);
fprintf('V_inner(U_match)   = %.8e\n', V_inner_match);
fprintf('jump inner - one   = %.8e\n', V_inner_match - V_one_match);

%%
% 7) Plotting parameters

col_two = [0.3010 0.7450 0.9330];
col_num = [0 0 0];
col_one = [1 0 0];

shift_num = -7.76073;
shift_two = -7.76073;
shift_one = -7.76073;

%%
% 8) Figure (b): U versus z
% Blue at bottom, red on top

fig1 = figure('Color','w','Position',[100 100 620 460]);
ax1 = axes(fig1);
hold(ax1,'on');

% Two-term outer part (plot first: bottom)
h_two = plot(ax1, z_two_outer + shift_two, U_two_outer, ...
    '-', 'Color', col_two, 'LineWidth', 1.6);

% Two-term inner part (plot first: bottom)
plot(ax1, z_two_inner + shift_two, U_two_inner, ...
    '-', 'Color', col_two, 'LineWidth', 2.2, ...
    'HandleVisibility', 'off');

% Numerical solution (middle layer)
h_num = plot(ax1, z_num + shift_num, U_num, ...
    '--', 'Color', col_num, 'LineWidth', 1.4);

% One-term solution (plot last: top)
h_one = plot(ax1, z_U0 + shift_one, U0_sol, ...
    '--', 'Color', col_one, 'LineWidth', 1.6);

xlabel(ax1, '$z$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax1, '$U$', 'Interpreter', 'latex', 'FontSize', 16);

legend(ax1, [h_two, h_num, h_one], ...
    {'Two-term solution', 'Numerical solution', 'One-term solution'}, ...
    'Interpreter', 'latex', 'Location', 'northeast', ...
    'FontSize', 9, 'Box', 'on');

xlim(ax1, [-16.5, 1]);
ylim(ax1, [-0.02, 1.02]);
grid(ax1, 'on');
box(ax1, 'on');
set(ax1, 'FontSize', 13, 'LineWidth', 1.0, ...
    'TickLabelInterpreter','latex');

%% Inset for travelling wave
ax1_inset = axes('Position', [0.24 0.28 0.25 0.25]);
hold(ax1_inset, 'on');

% Blue first
plot(ax1_inset, z_two_outer + shift_two, U_two_outer, ...
    '-', 'Color', col_two, 'LineWidth', 1.2);

plot(ax1_inset, z_two_inner + shift_two, U_two_inner, ...
    '-', 'Color', col_two, 'LineWidth', 1.8);

% Black second
plot(ax1_inset, z_num + shift_num, U_num, ...
    '--', 'Color', col_num, 'LineWidth', 1.1);

% Red last
plot(ax1_inset, z_U0 + shift_one, U0_sol, ...
    '--', 'Color', col_one, 'LineWidth', 1.2);

xlim(ax1_inset, [-16.2, 0.2]);
ylim(ax1_inset, [0, max(0.03, 1.2 * U_match)]);
grid(ax1_inset, 'on');
box(ax1_inset, 'on');
set(ax1_inset, 'FontSize', 8, 'LineWidth', 0.8, ...
    'TickLabelInterpreter','latex');

ax1_inset.XTick = [-16 -14 -12 -10 0];

%%
% 9) Figure (c): phase plane V versus U
% Blue at bottom, red on top

fig2 = figure('Color','w','Position',[180 120 620 460]);
ax2 = axes(fig2);
hold(ax2,'on');

% Two-term outer part (plot first: bottom)
h2_two = plot(ax2, U_two_outer_pp, V_two_outer_pp, ...
    '-', 'Color', col_two, 'LineWidth', 1.5);

% Two-term inner part (plot first: bottom)
plot(ax2, U_two_inner_pp, V_two_inner_pp, ...
    '-', 'Color', col_two, 'LineWidth', 2.4, ...
    'HandleVisibility', 'off');

% Numerical solution (middle)
h2_num = plot(ax2, U_num_pp, V_num_pp, ...
    '--', 'Color', col_num, 'LineWidth', 1.4);

% One-term solution (plot last: top)
h2_one = plot(ax2, U0_pp, V0_pp, ...
    '--', 'Color', col_one, 'LineWidth', 1.6);

xlabel(ax2, '$U$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax2, '$V$', 'Interpreter', 'latex', 'FontSize', 16);

legend(ax2, [h2_two, h2_num, h2_one], ...
    {'Two-term solution', 'Numerical solution', 'One-term solution'}, ...
    'Interpreter', 'latex', 'Location', 'northeast', ...
    'FontSize', 9, 'Box', 'on');

xlim(ax2, [0, 1]);
ylim(ax2, [-0.35, 0.01]);
grid(ax2, 'on');
box(ax2, 'on');
set(ax2, 'FontSize', 13, 'LineWidth', 1.0, ...
    'TickLabelInterpreter','latex');

%% Inset for small-U phase plane

ax2_inset = axes('Position', [0.18 0.26 0.30 0.30]);
hold(ax2_inset, 'on');

% Blue first
plot(ax2_inset, U_two_outer_pp, V_two_outer_pp, ...
    '-', 'Color', col_two, 'LineWidth', 1.1);

plot(ax2_inset, U_two_inner_pp, V_two_inner_pp, ...
    '-', 'Color', col_two, 'LineWidth', 2.2);

% Black second
plot(ax2_inset, U_num_pp, V_num_pp, ...
    '--', 'Color', col_num, 'LineWidth', 1.0);

% Red last
plot(ax2_inset, U0_pp, V0_pp, ...
    '--', 'Color', col_one, 'LineWidth', 1.2);

xlim(ax2_inset, [0, 1.2 * U_match]);

% automatic local ylim based on the plotted local data
local_mask_num = U_num_pp <= 1.2 * U_match;
local_mask_one = U0_pp <= 1.2 * U_match;
local_mask_in  = U_two_inner_pp <= 1.2 * U_match;

local_vals = [V_num_pp(local_mask_num); ...
              V0_pp(local_mask_one); ...
              V_two_inner_pp(local_mask_in)];

local_vals = local_vals(isfinite(local_vals));

if isempty(local_vals)
    ylim(ax2_inset, [-0.12, 0.01]);
else
    y_min_local = min(local_vals);
    y_max_local = max(local_vals);

    if y_min_local < 0
        ylim(ax2_inset, [1.15 * y_min_local, max(0.01, 1.15 * y_max_local)]);
    else
        ylim(ax2_inset, [-0.12, 0.01]);
    end
end

grid(ax2_inset, 'on');
box(ax2_inset, 'on');
set(ax2_inset, 'FontSize', 8, 'LineWidth', 0.8, ...
    'TickLabelInterpreter','latex');

% optional vector export
% exportgraphics(fig1,'fig9b_explicit_inner_blue_bottom.pdf','ContentType','vector');
% exportgraphics(fig2,'fig9c_explicit_inner_blue_bottom.pdf','ContentType','vector');

%%
% Local functions

function [z_sol, U0_sol, V0_sol] = solve_U0(F, c)
    % Solve one-term travelling wave ODE in z:
    % U' = V,  V' = -cV - F(U)

    odefun = @(z, Y) [Y(2); -c*Y(2) - F(Y(1))];

    z_span = linspace(-12, 50, 2500);
    Y0 = [1; -1e-4];

    opts = odeset( ...
        'RelTol', 1e-9, ...
        'AbsTol', 1e-12, ...
        'MaxStep', 5e-3);

    [z_sol, Y_sol] = ode15s(odefun, z_span, Y0, opts);

    U0_sol = Y_sol(:,1);
    V0_sol = Y_sol(:,2);

    keep = isfinite(U0_sol) & isfinite(V0_sol) & ...
           (U0_sol >= -0.05) & (U0_sol <= 1.05);

    z_sol  = z_sol(keep);
    U0_sol = U0_sol(keep);
    V0_sol = V0_sol(keep);
end

function [U_sol, V_combined] = solve_combined_V(F, c, epsilon)
    % Solve phase-plane equation:
    % dV/dU = (-cV - F(U))/V
    %
    % and first-order correction V = V0 + epsilon V1

    U_values = linspace(1, 1e-6, 2500).';

    dV0dU = @(U, V0) (-c*V0 - F(U)) ./ V0;

    V0_start = -epsilon;
    opts = odeset( ...
        'RelTol', 1e-9, ...
        'AbsTol', 1e-12, ...
        'MaxStep', 5e-3);

    [U_sol_V0, V0_sol] = ode15s(dV0dU, U_values, V0_start, opts);

    U_sol_V0 = U_sol_V0(:);
    V0_sol   = V0_sol(:);

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

    V1_interp = interp1(U_sol_V1, V1_sol, U_sol_V0, ...
        'linear', 'extrap');

    V_combined = V0_sol + epsilon * V1_interp;
    U_sol = U_sol_V0;

    keep = isfinite(U_sol) & isfinite(V_combined);
    U_sol = U_sol(keep);
    V_combined = V_combined(keep);
end

function yq = interp_monotone(x, y, xq)
    % Safe interpolation for either increasing or decreasing x.

    x = x(:);
    y = y(:);

    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);

    if isempty(x)
        error('interp_monotone: empty x array.');
    end

    % Sort in ascending x and remove duplicate x-values.
    [x_sort, idx] = sort(x, 'ascend');
    y_sort = y(idx);

    [x_unique, ia] = unique(x_sort, 'stable');
    y_unique = y_sort(ia);

    yq = interp1(x_unique, y_unique, xq, 'linear', 'extrap');
end
clear; clc; close all;

%% =========================
% Parameters
% =========================
a = 0.3;                                 % Allee threshold
c_star = sqrt(2) * (0.5 - a);            % critical speed
epsilon = 1e-3;                          % small parameter
c = c_star - epsilon;

F = @(U) U .* (1 - U) .* (U - a);

%% =========================
% Matching point
% Use a larger U_match first so that the inner branch is visible.
% If needed, change this back to:
% U_match = 5 * epsilon;
%% =========================
U_match = 0.02;
% U_match = 5 * epsilon;

err_floor = 1e-14;   % for semilogy stability

%% =========================
% 1) One-term solution V0(U)
% =========================
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

%% =========================
% 2) Numerical solution V_num(U)
% =========================
[U_num, V_num] = solve_combined_V(F, c_star, epsilon);

U_num = U_num(:);
V_num = V_num(:);

% Ensure U decreases from 1 to 0
if U_num(1) < U_num(end)
    U_num = flipud(U_num);
    V_num = flipud(V_num);
end

%% =========================
% 3) Explicit inner small-U branch
%
% psi_0(phi) ~ -gamma * phi + D * phi^(-r)
%
% U = epsilon^p * phi
% V = epsilon^q * psi
%
% where
% p = q = r = 1/(2a+1),
% gamma = 1/sqrt(2),
% D = -Gamma(2a+1) Gamma(3-2a) / 6
%% =========================
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
valid_idx = isfinite(U_small) & isfinite(V_small) & ...
            (U_small > 0) & (V_small < -1e-12);

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

%% =========================
% 4) Build two-term phase-plane curve V_two(U)
%
% U >= U_match : one-term solution
% U <  U_match : explicit inner branch
%% =========================

% One-term value at matching point
V_one_match = interp_monotone(U0_sol, V0_sol, U_match);

% Inner value at matching point
V_inner_match = interp_monotone(U_small, V_small, U_match);

% ---------- outer part: one-term branch ----------
idx_outer = U0_sol > U_match;

U_two_outer = [U0_sol(idx_outer); U_match];
V_two_outer = [V0_sol(idx_outer); V_one_match];

% ---------- inner part: explicit inner asymptotic branch ----------
idx_inner = U_small < U_match;

U_two_inner = [U_match; U_small(idx_inner)];
V_two_inner = [V_inner_match; V_small(idx_inner)];

% ---------- combine outer + inner ----------
% remove duplicate U_match from the inner part
U_two = [U_two_outer; U_two_inner(2:end)];
V_two = [V_two_outer; V_two_inner(2:end)];

%% =========================
% 5) Sort all curves as ascending U for interpolation
% =========================
[U_num_pp, idx_num_pp] = sort(U_num, 'ascend');
V_num_pp = V_num(idx_num_pp);

[U_one_pp, idx_one_pp] = sort(U0_sol, 'ascend');
V_one_pp = V0_sol(idx_one_pp);

[U_two_pp, idx_two_pp] = sort(U_two, 'ascend');
V_two_pp = V_two(idx_two_pp);

%% =========================
% 6) Build common U-grid for error comparison
%
% We compare on the common overlap domain of all three curves.
%% =========================
U_min_common = max([min(U_num_pp), min(U_one_pp), min(U_two_pp)]);
U_max_common = min([max(U_num_pp), max(U_one_pp), max(U_two_pp)]);

if U_min_common >= U_max_common
    error('No common U interval exists for the error comparison.');
end

N_eval = 4000;
U_eval = linspace(U_min_common, U_max_common, N_eval).';

V_num_eval = interp_monotone(U_num_pp, V_num_pp, U_eval);
V_one_eval = interp_monotone(U_one_pp, V_one_pp, U_eval);
V_two_eval = interp_monotone(U_two_pp, V_two_pp, U_eval);

%% =========================
% 7) Errors
% =========================
err_one = abs(V_num_eval - V_one_eval);
err_two = abs(V_num_eval - V_two_eval);

err_one_log = max(err_one, err_floor);
err_two_log = max(err_two, err_floor);

fprintf('\nDiagnostics:\n');
fprintf('Common U range     = [%.6e, %.6e]\n', U_min_common, U_max_common);
fprintf('Number eval points = %d\n', N_eval);
fprintf('max err (one-term) = %.8e\n', max(err_one));
fprintf('max err (two-term) = %.8e\n', max(err_two));
fprintf('min err (one-term) = %.8e\n', min(err_one));
fprintf('min err (two-term) = %.8e\n', min(err_two));

%% =========================
% 8) Figure 1: Error plot (linear scale)
%% =========================
fig1 = figure('Color','w','Position',[100 100 620 460]);
ax1 = axes(fig1);
hold(ax1,'on');

h_two = plot(ax1, U_eval, err_two, ...
    '-', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.8);

h_one = plot(ax1, U_eval, err_one, ...
    'r--', 'LineWidth', 1.6);

xlabel(ax1, '$U$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax1, '$|V_{\mathrm{num}}(U)-V_{\mathrm{asym}}(U)|$', ...
    'Interpreter', 'latex', 'FontSize', 16);

legend(ax1, [h_two, h_one], ...
    {'Two-term error', 'One-term error'}, ...
    'Interpreter', 'latex', 'Location', 'northeast', ...
    'FontSize', 9, 'Box', 'on');

xlim(ax1, [U_min_common, U_max_common]);
grid(ax1, 'on');
box(ax1, 'on');
set(ax1, 'FontSize', 13, 'LineWidth', 1.0, ...
    'TickLabelInterpreter','latex');

%% =========================
% 9) Figure 2: Error plot (semilogy)
%% =========================
fig2 = figure('Color','w','Position',[180 120 620 460]);
ax2 = axes(fig2);
hold(ax2,'on');

h2_two = semilogy(ax2, U_eval, err_two_log, ...
    '-', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.8);

h2_one = semilogy(ax2, U_eval, err_one_log, ...
    'r--', 'LineWidth', 1.6);

xlabel(ax2, '$U$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax2, '$|V_{\mathrm{num}}(U)-V_{\mathrm{asym}}(U)|$', ...
    'Interpreter', 'latex', 'FontSize', 16);

legend(ax2, [h2_two, h2_one], ...
    {'Two-term error', 'One-term error'}, ...
    'Interpreter', 'latex', 'Location', 'northeast', ...
    'FontSize', 9, 'Box', 'on');

xlim(ax2, [U_min_common, U_max_common]);
grid(ax2, 'on');
box(ax2, 'on');
set(ax2, 'FontSize', 13, 'LineWidth', 1.0, ...
    'TickLabelInterpreter','latex');

% optional export
% exportgraphics(fig1,'error_plot_linear.pdf','ContentType','vector');
% exportgraphics(fig2,'error_plot_semilogy.pdf','ContentType','vector');

%% =========================
% Local functions
% =========================
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
    % Safe interpolation for monotone data

    x = x(:);
    y = y(:);

    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);

    if isempty(x)
        error('interp_monotone: empty x array.');
    end

    % Sort ascending and remove duplicates
    [x_sort, idx] = sort(x, 'ascend');
    y_sort = y(idx);

    [x_unique, ia] = unique(x_sort, 'stable');
    y_unique = y_sort(ia);

    yq = interp1(x_unique, y_unique, xq, 'linear', 'extrap');
end
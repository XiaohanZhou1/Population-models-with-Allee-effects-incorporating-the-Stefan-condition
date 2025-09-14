% Parameters
a = 0.3;
c = 1e-3;
F = @(U) U .* (1 - U) .* (U - a);  % Reaction term

% Numerical solution by ODE
ode = @(U, V) -c - F(U)./V;

% Integrate from U=1 to U=0
epsilon = 1e-8; % Small perturbation to start integration
U_range = [1, 0]; % Integration interval
V_init = -epsilon; % V(1) ~ 0

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[U_num, V_num] = ode15s(ode, U_range, V_init, options);

% One-term asymptotic solution
V0_integral = @(U) arrayfun(@(u) integral(@(s) F(s), u, 1, 'RelTol',1e-8,'AbsTol',1e-10), U);
V0 = @(U) -sqrt(2 * V0_integral(U));

U_vals = linspace(0, 1, 500);  % Uniform domain for plotting
V0_vals = V0(U_vals);

% Two-term asymptotic solution
% Define V1(U) = [∫_U^1 V0(s) ds] / V0(U)
V1_integral = @(U) arrayfun(@(u) integral(@(s) V0(s), u, 1, 'RelTol',1e-8,'AbsTol',1e-10), U);
V1 = @(U) V1_integral(U) ./ V0(U);
V_two_term = @(U) V0(U) + c * V1(U);
V2_vals = V_two_term(U_vals);

% Interpolate numerical to match U_vals for comparison
V_num_interp = interp1(U_num, V_num, U_vals, 'linear', 'extrap');

figure;
plot(U_vals, V_num_interp, 'b-', 'LineWidth', 2, 'DisplayName', 'Numerical (ODE)');
hold on;
plot(U_vals, V0_vals, 'r-', 'LineWidth', 2, 'DisplayName', 'One-term Asymptotic');
plot(U_vals, V2_vals, 'g--', 'LineWidth', 2, 'DisplayName', 'Two-term Asymptotic');
xlabel('U');
ylabel('V');
legend('Location', 'Best');
grid on;
hold off;

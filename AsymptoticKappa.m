% Numerical compare
%kappa epsilon
clc; clear;

% Asymptotic
a = 0.25;
c_star = sqrt(2) * (0.5 - a);
gamma_val = -1 / sqrt(2);
rho = 1 / (sqrt(2) * (2*a + 1));
p = 1 / (2*a + 1);

D = -gamma(2*a + 1) * gamma(3 - 2*a) / gamma(4);
xi = (- gamma_val + rho)^rho * (-D)^(-gamma_val);
y0_0 = -xi^(1 / (rho - gamma_val));

eps1 = logspace(-8, 0, 100);
kappa_match = -(c_star) ./ y0_0 .* eps1.^(-p);

% Full phase plane
eps2 = logspace(-8, 0, 100);
kappa_wave = zeros(size(eps2));
V_at_0_list = zeros(size(eps2));
F = @(U) U .* (1 - U) .* (U - a);
z_span = linspace(-55, 0, 2000);
delta = 1e-9;
options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.01);

for i = 1:length(eps2)
    epsilon = eps2(i);
    c = c_star - epsilon;

    odefun = @(z, Y) [
        Y(2); 
        -(c * Y(2) + F(Y(1)))
    ];

    Y0 = [1 - delta; -delta];
    [~, Y_sol] = ode15s(odefun, z_span, Y0, options);
    [~, idx_zero] = min(abs(Y_sol(:,1)));
    V_at_0 = Y_sol(idx_zero, 2);
    V_at_0_list(i) = V_at_0;
    kappa_wave(i) = - (c_star - epsilon) / V_at_0;
end

% Plot
figure;

loglog(eps2, kappa_wave, 'b--', 'LineWidth', 2); hold on;
loglog(eps1, kappa_match, 'r-', 'LineWidth', 2); 
xlabel('$\epsilon$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', 20);
legend({'Numerical TW', 'Large-$c$ expansion'}, 'Location', 'best', ...
    'Interpreter', 'latex');
grid on;




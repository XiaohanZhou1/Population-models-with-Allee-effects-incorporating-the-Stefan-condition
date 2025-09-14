a = 0.25;
c_star = sqrt(2) * (0.5 - a);
gamma_val = -1 / sqrt(2); % gamma < 0
rho = 1 / (sqrt(2) * (2*a + 1)); % rho > 0
epsilon = 1e-8;
c = c_star - epsilon;

% gamma function based D/E
D = -gamma(2*a + 1) * gamma(3 - 2*a) / gamma(4);
E = D;

% p, xi, y0(0), kappa
p = 1 / (2*a + 1);
xi = (- gamma_val + rho)^rho * (-D)^(-gamma_val);
y0_0 = -xi^(1 / (rho - gamma_val));

implicit_fun = @(y0, x) abs(y0 - rho*x)^rho * abs(y0 - gamma_val*x)^(-gamma_val) - xi;

x_vals = linspace(0, 10, 500);
y_vals = NaN(size(x_vals));
y_vals(1) = y0_0;

for i = 2:length(x_vals)
    x = x_vals(i);
    y_prev = y_vals(i-1);

    fun = @(y0) implicit_fun(y0, x);

    if isnan(y_prev)
        y_guess = gamma_val * x - 0.01;
        try
            y_vals(i) = fzero(fun, y_guess);
        catch
            y_vals(i) = NaN;
        end
    else
        try
            y_low = y_prev - 0.1;
            y_high = y_prev + 0.1;
            y_vals(i) = fzero(fun, [y_low, y_high]);
        catch
            y_vals(i) = NaN;
        end
    end
end

figure;
plot(x_vals, y_vals, 'b-', 'LineWidth', 2); hold on;
plot(x_vals, gamma_val * x_vals, 'g--', 'LineWidth', 1.2); % y = gamma x
plot(0, y0_0, 'ko', 'MarkerFaceColor', 'k', 'HandleVisibility','off');

xlabel('x');
ylabel('y_0');
legend('y_0(x)', 'y = \gamma x', 'y_0(0)', 'Location', 'northwest');
xlim([0, 2.5]);
ylim([-2, 0]);
grid on;

% y0(x) = gamma*x + D*x^(rho/gamma)
y_extra = gamma_val * x_vals + D * x_vals.^(rho / gamma_val);

plot(x_vals, y_extra, 'r-.', 'LineWidth', 1.5);

legend('y_0(x)', 'y = \gamma x', 'y_0(x) = \gamma x + D x^{\rho/\gamma}', ...
       'Location', 'northwest');


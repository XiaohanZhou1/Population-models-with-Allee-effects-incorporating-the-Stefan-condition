% Numerically simulated wave speeds
kappa_values = [20, 30, 40];
a_values = linspace(0, 0.5, 8);

% Store the results
results = zeros(length(kappa_values), length(a_values));

for i = 1:length(kappa_values)
    for j = 1:length(a_values)
        c = AlleeStefanPDE(kappa_values(i), a_values(j));
        results(i, j) = c;  % Assume this function returns a speed c
    end
end

% Theoretically simulated wave speeds
a_theory = linspace(0, 0.55, 10000);
c_theory = sqrt(2).*(1/2-a_theory);  % Simplified theoretical curve

% Prepare the figure
figure;
colors = lines(length(kappa_values));

% Plot numerical results
hold on;
for i = 1:length(kappa_values)
    % Display Name uses LaTeX format for kappa
    plot(a_values, results(i, :), 'o-', 'DisplayName', sprintf('$\\kappa = %d$', kappa_values(i)), 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', colors(i, :));
end

% Plot theoretical results
plot(a_theory, c_theory, 'k-', 'LineWidth', 2, 'DisplayName', 'Theoretical');

% Set the axis labels and legend
xlabel('a');
ylabel('c');
legend show;
set(legend, 'Interpreter', 'latex');  % Ensure the legend uses LaTeX interpreter
grid on;
box on;
xlim([0, 0.5]);
ylim([0, 0.8]);

hold off;



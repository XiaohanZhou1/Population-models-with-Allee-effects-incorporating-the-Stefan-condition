cValues=[0.5, 1, 2, 10];
k=1;
    
figure;
hold on;

options = odeset('RelTol', 1e-12);

for c = cValues
    ode = @(z, Y) [Y(2); -c*Y(2)-Y(1)*(1-Y(1))];
        
    z = [-1000, 0];
    Y0 = [1; -0.001]; % Initial point close to (1, 0)

    [Z, Y] = ode15s(ode, z, Y0, options);

    plot(Y(:,1), Y(:,2), 'DisplayName', ['c = ', num2str(c)], 'LineWidth', 1);
end

xlabel('U');
ylabel('V');
legend show;
grid on;
box on;
hold off;

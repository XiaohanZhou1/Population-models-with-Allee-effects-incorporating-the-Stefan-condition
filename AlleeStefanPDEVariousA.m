function AlleeStefanPDEVariousA(kappa)
    Y = 1; % Y domain corresponds to the transformed X domain
    N1 = 2e+3;
    dy = Y/(N1-1);
    y = linspace(0, Y, N1);

    T = 100;
    N2 = 1e+4;
    dt = T/N2;

    M = 1e+5; % M should be sufficiently large
    L0 = 1; % Initial value of L(t) at t=0
    L = L0; % Initialize L(t) with its initial value
    alpha = M + L; % Initial value of M at t=0

    p.dy = dy;
    p.dt = dt;
    p.y = y(:);
    p.L = L;
    p.alpha = alpha;
    p.M = M;
    p.kappa = kappa;
    
    a_values = [0.3, 0.5, 0.7];  % Array of a values
    colors = ['r', 'g', 'b'];    % Colors for each plot
    lineStyles = {'-', '--', ':'};  % Line styles
    
    figure; hold on;  % Prepare figure to hold all plots
    
    for idx = 1:length(a_values)
        p.a = a_values(idx);  % Set current a value
        u = zeros(N1, 1);
        u(:, 1) = 1 * (1 - (y > 1));  % Reset initial condition
        
        % RK4 time-stepping for current a value
        for n = 1:N2
            du_dy_at_y1 = (u(N1) - u(N1 - 1)) / p.dy;
            p.c = -p.kappa / p.alpha * du_dy_at_y1;
            
            u = RK4total(p, u);
            
            % Update L using the Stefan condition
            du_dy_at_L = (u(N1) - u(N1 - 1)) / p.dy;
            p.L = p.L - (p.kappa / p.alpha) * du_dy_at_L * p.dt;
    
            p.alpha = p.M + p.L;
    
            u(1) = u(2); % Apply Neumann boundary condition du/dy=0 at y=0
            u(end) = 0; % Apply Dirichlet boundary condition u(y=1,t)=0
            
            % Debugging output to check the value of a
            if mod(n, 1000) == 0
                disp(['Step: ', num2str(n), ', a = ', num2str(p.a)]);
            end
        end
    
        x = -p.M + p.y .* (p.L + p.M);

        U = u;
        z = x - p.c * T;
        I = find(U < 1e-6, 1, 'first');
        z0 = z(I);
    
        plot(x - z0, U, 'LineWidth', 2, 'Color', colors(idx), 'LineStyle', lineStyles{idx}); % Plot with specific color and line style
        pause(0.01); % Pause to update
    end
    
    xlim([-30, 30]);
    ylim([0, 1]);
    xlabel('z');
    ylabel('U');
    legend(arrayfun(@(a) ['a = ' num2str(a)], a_values, 'UniformOutput', false));
    box on;
    hold off;
end

% Remaining functions (RK4total and f_RK4) remain unchanged, except for using p.a instead of a hardcoded a value
function [dx] = RK4total(p, u)
    k1 = p.dt * f_RK4(p, 0, u);
    k2 = p.dt * f_RK4(p, k1 / 2, u);
    k3 = p.dt * f_RK4(p, k2 / 2, u);
    k4 = p.dt * f_RK4(p, k3, u);
    dx = u + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
end

function f = f_RK4(p, ~, u)
    dy = p.dy;
    alpha = p.alpha;

    du_dy = [0; (u(3:end) - u(1:end-2)) / (2 * dy); 0];
    d2u_dy2 = [0; (u(3:end) - 2 * u(2:end-1) + u(1:end-2)) / (dy^2); 0];
    
    f = alpha^(-1) * du_dy .* (p.c .* p.y) + alpha^(-2) * d2u_dy2 + u .* (1 - u) .* (u - p.a); % RHS of the PDE
end

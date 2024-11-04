function compareWaves(kappa, a) % Figure 9(c)
    % First run
    u_final_1 = solvePDE(kappa, a, []);
    
    % Second run with the final state of the first run as initial condition
    u_final_2 = solvePDE(kappa, a, u_final_1);
    
    % Plot the traveling waves for comparison
    plotWaves(u_final_1, u_final_2);
    print('-depsc2', 'figure.eps');
end

function u_final = solvePDE(kappa, a, initial_u)
    Y = 1; % Y domain corresponds to the transformed X domain
    N1 = 2e+5;
    dy = Y/(N1-1);
    y = linspace(0, Y, N1);

    T = 100;
    N2 = 1e+4;
    dt = T/N2;

    M = 1e+5; % M should be sufficiently large
    L0 = 1; % Initial value of L(t) at t=0
    L = L0; % Initialize L(t) with its initial value
    alpha = M+L; % Initial value of M at t=0

    u = zeros(N1, 1); % The matrix that includes numerical solutions U(y, t)
    if isempty(initial_u)
        u(:, 1) = 1*(1-(y>1)); % Setting initial condition based on the Heaviside function
    else
        u(:, 1) = initial_u; % Use the provided initial condition
    end
    
    p.dy = dy;
    p.dt = dt;
    p.y = y(:);
    p.L = L;
    p.alpha = alpha;
    p.M = M;
    p.kappa = kappa;
    p.a = a;
    
    % RK4 time-stepping
    for n = 1:N2
        du_dy_at_y1 = (u(N1)-u(N1-1))/p.dy;
        p.c = -p.kappa/p.alpha*du_dy_at_y1;
        
        u = RK4total(p, u);
        
        % Update L using the Stefan condition
        du_dy_at_L = (u(N1)-u(N1-1))/p.dy;
        p.L = p.L - (p.kappa/p.alpha) * du_dy_at_L * p.dt;

        p.alpha = p.M + p.L;

        u(1) = u(2); % Apply Neumann boundary condition du/dy=0 at y=0
        u(end) = 0; % Apply Dirichlet boundary condition u(y=1,t)=0
    end

    u_final = u(:, end); % Return the final state of u

    x = -p.M + p.y .* (p.L + p.M);
    U = u_final;
    z = x - p.c * T;
    I = find(U < 1e-6, 1, 'first');
    z0 = z(I);
end

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

    % Forward finite difference with second order accuracy
    du_dy = zeros(size(u));
    du_dy(1:end-2) = (-3 * u(1:end-2) + 4 * u(2:end-1) - u(3:end)) / (2 * dy);
    du_dy(end-1:end) = (u(end-1:end) - u(end-2:end-1)) / dy;

    d2u_dy2 = [0; (u(3:end) - 2 * u(2:end-1) + u(1:end-2)) / (dy^2); 0];
    
    f = alpha^(-1) * du_dy .* (p.c .* p.y) + alpha^(-2) * d2u_dy2 + u .* (1 - u) .* (u - p.a); % RHS of the PDE
end

function plotWaves(u1, u2)
    y_shift = linspace(0, 2e+5, 2e+5) - 2e+5; % y axis shifted by -2e+5

    figure;
    hold on;
    plot(y_shift, u1, 'LineWidth', 2);
    plot(y_shift, u2, '--', 'LineWidth', 2); % Use dashed line for the second run
    legend('First Run', 'Second Run');
    xlabel('y');
    ylabel('u');
    xlim([-25, 0]);
    ylim([0, 1]);
    box on;
    grid off;
    hold off;
end

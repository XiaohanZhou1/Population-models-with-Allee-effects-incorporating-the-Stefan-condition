function c = AlleeStefanPDESurfacePlot(kappa, a)%y-t
    % Parameters
    Y = 1; % Y domain corresponds to the transformed X domain
    N1 = 2e+3; % Initial size for debugging
    dy = Y / (N1 - 1);
    y = linspace(0, Y, N1);

    T = 10; % Adjust as needed
    N2 = 1e+3; % Initial size for debugging
    dt = T / N2;

    M = 1e+3; % Initial size for debugging
    L0 = 1; % Initial value of L(t) at t=0
    L = L0; % Initialize L(t) with its initial value
    alpha = M + L; % Initial value of alpha at t=0

    % Initial conditions
    u = zeros(N1, N2); % The matrix that includes numerical solutions U(y, t)
    u(:, 1) = 1 * (1 - (y > 1)); % Setting initial condition based on the Heaviside function
    
    % Parameters structure
    p.dy = dy;
    p.dt = dt;
    p.y = y(:);
    p.L = L;
    p.alpha = alpha;
    p.M = M;
    p.kappa = kappa;
    p.a = a;

    % Pre-allocate memory for wavefront position
    L_pos = zeros(1, N2);

    % RK4 time-stepping
    for n = 1:N2-1
        du_dy_at_y1 = (u(N1, n) - u(N1-1, n)) / p.dy;
        p.c = -p.kappa / p.alpha * du_dy_at_y1;
        
        u(:, n+1) = RK4total(p, u(:, n));
        
        p.L = p.L + p.c * p.dt;
        p.alpha = p.M + p.L;

        % Record wavefront position
        L_pos(n+1) = p.L;

        % Debugging output
        if mod(n, 100) == 0
            disp(['Time step: ', num2str(n), ', L(t): ', num2str(p.L), ', c: ', num2str(p.c)]);
        end

        u(1, n+1) = u(2, n+1); % Apply Neumann boundary condition du/dy=0 at y=0
        u(end, n+1) = 0; % Apply Dirichlet boundary condition u(y=1,t)=0
    end

    % Compute wave speed as the change in wavefront position over time
    t = linspace(0, T, N2);
    valid_indices = L_pos > 0; % Only consider valid L positions
    p_fit = polyfit(t(valid_indices), L_pos(valid_indices), 1);
    c = p_fit(1); % The slope of the linear fit represents the wave speed
    disp(['Estimated wave speed c = ', num2str(c)]);
    
    % Convert to double precision for surf
    y = double(y);
    t = double(t);
    u = double(u);

    % Debugging sizes
    disp(['Size of y: ', num2str(size(y))]);
    disp(['Size of t: ', num2str(size(t))]);
    disp(['Size of u: ', num2str(size(u))]);

    % Plot the 2D surface plot
    figure;
    surf(t, y, u, 'EdgeColor', 'none');
    view(2); % View from above
    colorbar;
    box on;
    grid off;
    xlabel('t');
    ylabel('y');
    
    % Plot the 3D surface plot (y, t, u)
    figure;
    [T, Y] = meshgrid(t, y);
    surf(T, Y, u, 'EdgeColor', 'none'); % Note: u instead of u'
    set(gca, 'YDir', 'reverse'); % Reverse the direction of the y-axis
    colorbar;
    box on;
    grid off;
    xlabel('t');
    ylabel('y');
    zlabel('u');
    
    % Additional debugging: display the range of u
    disp(['Range of u: ', num2str(min(u(:))), ' to ', num2str(max(u(:)))]);
end

function u_next = RK4total(p, u)
    k1 = p.dt * f_RK4(p, u);
    k2 = p.dt * f_RK4(p, u + 0.5 * k1);
    k3 = p.dt * f_RK4(p, u + 0.5 * k2);
    k4 = p.dt * f_RK4(p, u + k3);
    u_next = u + (k1 + 2 * k2 + 2 * k3 + k4) / 6;  % RK4 step
end

function f = f_RK4(p, u)
    dy = p.dy;
    alpha = p.alpha;
    a = p.a;

    % Computing derivatives
    du_dy = zeros(size(u));
    d2u_dy2 = zeros(size(u));

    du_dy(2:end-1) = (u(3:end) - u(1:end-2)) / (2 * dy);
    d2u_dy2(2:end-1) = (u(3:end) - 2 * u(2:end-1) + u(1:end-2)) / (dy^2);

    % Ensure p.y is a column vector
    if size(p.y, 2) > 1
        p.y = p.y';
    end

    % Calculating the function f for RK4
    f = alpha^(-1) * du_dy .* (p.c * p.y) + alpha^(-2) * d2u_dy2 + u .* (1 - u) .* (u - a); % RHS of the PDE
end







function FisherStefanPlotErrors
    kappa_values = linspace(10, 30, 20);
    kappa_calculated_values = zeros(size(kappa_values));
    relative_errors = zeros(size(kappa_values));

    for i = 1:length(kappa_values)
        % Input kappa, solve the PDE to obtain c
        kappa = kappa_values(i);

        Y = 1; % Y domain corresponds to the transformed X domain
        N1 = 1e+5; % Increased size for better accuracy
        dy = Y/(N1-1);
        y = linspace(0, Y, N1);

        T = 10;
        N2 = 1e+4; % Increased size for better accuracy
        dt = T/N2;

        M = 10000; % M should be sufficiently large
        L0 = 1; % Initial value of L(t) at t=0
        L = L0; % Initialize L(t) with its initial value
        alpha = M+L; % Initial value of M at t=0

        u = zeros(N1, 1); % The matrix that includes numerical solutions U(y, t)
        u(:, 1) = 0.5*(1-(y>1)); % Setting initial condition based on the Heaviside function

        p.dy = dy;
        p.dt = dt;
        p.y = y(:);
        p.L = L;
        p.alpha = alpha;
        p.M = M;
        p.kappa = kappa;

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

        c = p.c;
        disp(['Estimated wave speed c = ', num2str(c)]);

        % Input c, solve the ODE system to obtain kappa
        try
            kappa_calculated = PhasePlane(c);
        catch ME
            disp(['ODE solver failed: ', ME.message]);
            kappa_calculated = NaN;
        end

        kappa_calculated_values(i) = kappa_calculated;
        relative_errors(i) = abs(kappa - kappa_calculated)/kappa;
    end

    % Generate the plot
    plot(kappa_values, relative_errors, '-o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'LineWidth', 2);
    xlabel('\kappa');
    ylabel('Relative Error');
    grid on;
    box on;
end

function dydz = odes(z, y, c)
    U = y(1);
    V = y(2);
    dUdz = V;
    dVdz = -c*V-U*(1-U);
    dydz = [dUdz; dVdz];
end

function [dx] = RK4total(p, u)
    k1 = p.dt*f_RK4(p,0,u);
    k2 = p.dt*f_RK4(p,k1/2,u);
    k3 = p.dt*f_RK4(p,k2/2,u);
    k4 = p.dt*f_RK4(p,k3,u);
    dx = u + (k1+2*k2+2*k3+k4)/6;
end

function f = f_RK4(p, ~, u)
    dy = p.dy;
    alpha = p.alpha;

    % Forward finite difference with second order accuracy
    du_dy = zeros(size(u));
    du_dy(1:end-2) = (-3*u(1:end-2) + 4*u(2:end-1) - u(3:end)) / (2*dy);
    du_dy(end-1:end) = (u(end-1:end) - u(end-2:end-1)) / dy;

    d2u_dy2 = [0; (u(3:end) - 2*u(2:end-1) + u(1:end-2))/(dy^2); 0];

    f = alpha^(-1)*du_dy.*(p.c.*p.y)+alpha^(-2)*d2u_dy2+u.*(1-u); % RHS of the PDE
end

function kappa = PhasePlane(c)
    ode = @(U, V) -c-U*(1-U)/V;

    epsilon = 1e-8;
    u = [1, 0]; % We integrate from U=1 to U=0
    V1 = -epsilon; % When U=1, V=0
    
    options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'MaxStep', 1e-3); % Increased accuracy
    [U, V] = ode15s(ode, u, V1, options); % Solve the ODE
    
    % Check if there are sufficient data points for interpolation
    if length(U) < 2
        error('Insufficient data points for interpolation');
    end

    V0 = interp1(U, V(:,1), 0, 'linear', 'extrap'); % The value of V at U=0
    kappa = -c/V0;
    
    % disp(['Calculated kappa: ', num2str(kappa)]);
end




  

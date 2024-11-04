function c = AlleeStefanPlotLt(kappa, a)
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
    u(:, 1) = 1*(1-(y>1)); % Setting initial condition based on the Heaviside function
    
    p.dy = dy;
    p.dt = dt;
    p.y = y(:);
    p.L = L;
    p.alpha = alpha;
    p.M = M;
    p.kappa = kappa;
    p.a = a;

    L_values = zeros(N2, 1);
    L_values(1) = L0;
    
    % RK4 time-stepping
    for n = 1:N2
        du_dy_at_y1 = (u(N1)-u(N1-1))/p.dy;
        p.c = -p.kappa/p.alpha*du_dy_at_y1;
        
        u = RK4total(p, u);
        
        % Update L using the Stefan condition
        du_dy_at_L = (u(N1)-u(N1-1))/p.dy;
        p.L = p.L - (p.kappa/p.alpha) * du_dy_at_L * p.dt;

        L_values(n) = p.L;

        p.alpha = p.M+p.L;

        u(1) = u(2); % Apply Neumann boundary condition du/dy=0 at y=0
        u(end) = 0; % Apply Dirichlet boundary condition u(y=1,t)=0
    end

    c = p.c;
    disp(['Estimated wave speed c = ', num2str(c)]);

    t = linspace(0, T, N2);
    plot(t, L_values, 'LineWidth', 2);
    hold on;
    
    % Second order derivative calculation using central difference
    L_prime_second_order = zeros(N2-2, 1);
    for n = 2:N2-1
        L_prime_second_order(n-1) = (L_values(n+1) - L_values(n-1)) / (2*dt);
    end
    
    % Find the linear region where the slope doesn't change much
    tolerance = 1e-6; % Tolerance for slope change
    linear_index = find(abs(diff(L_prime_second_order)) < tolerance, 1);
    
    if ~isempty(linear_index)
        linear_slope = L_prime_second_order(linear_index);
        disp(['Slope in the linear region = ', num2str(linear_slope)]);
        
        % Find the point in time t where the linear region starts
        t_linear_start = t(linear_index);
        
        % Plot the linear extension line
        L_linear_ext = linear_slope * (t - t_linear_start) + L_values(linear_index);
        plot(t, L_linear_ext, '--', 'LineWidth', 2);
        
        % Calculate the intersection with the y-axis
        intersection_y = L_values(linear_index) - linear_slope * t_linear_start;
        disp(['Intersection with y-axis = ', num2str(intersection_y)]);
    else
        disp('No linear region found with the given tolerance.');
    end
    
    xlabel('t');
    ylabel('L');
    xlim([0, 20]);
    ylim([0, 15]);
    box on;
    hold off;

    print('-depsc2', 'figure.eps');
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
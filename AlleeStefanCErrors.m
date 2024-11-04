c_values = linspace(0.1, 0.5, 10);
kappa_calculated_values = zeros(size(c_values));
relative_errors = zeros(size(c_values));

for i = 1:length(c_values) 
    % Input c, solve the ODE system to obtain kappa
    c = c_values(i);
    
    z1 = -1e+13;
    z2 = 0;
    
    epsilon = 1e-10;
    U0 = 1-epsilon;
    V0 = -epsilon;
    y0 = [U0; V0];
    
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [Z, Y] = ode15s(@(z, y) odes(z, y, c), [z1, z2], y0, options);


    % Extract U and V from the solution Y
    U = Y(:,1);
    V = Y(:,2);
    
    % Find the index I where U is closest to zero
    [~, I] = min(abs(U));
    
    % Check if the corresponding V is close enough to be considered as V at U=0
    if abs(U(I)) < 1e-1
        % Calculate kappa using the formula kappa = -c/V0
        V0_at_U0 = V(I);  % V value when U is closest to zero
        kappa = -c / V0_at_U0;
        disp(['Calculated kappa = ', num2str(kappa)]);
    else
        disp('U=0 was not found in the calculated solution.');
    end


    kappa_calculated_values(i) = kappa;
    

    
    Y = 1; % Y domain corresponds to the transformed X domain
    N1 = 2e+3;
    dy = Y/(N1-1);
    y = linspace(0, Y, N1);
    
    T = 100;
    N2 = 1e+4;
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
    
        p.L = p.L+p.c*p.dt;
        p.alpha = p.M+p.L;
    
        u(1) = u(2); % Apply Neumann boundary condition du/dy=0 at y=0
        u(end) = 0; % Apply Dirichlet boundary condition u(y=1,t)=0
    end
    
    c_actual = p.c;
    disp(['Estimated wave speed c = ', num2str(c_actual)]);
    relative_errors(i) = abs(c - c_actual) / c;
end


% Generate the plot for relative errors
plot(c_values, relative_errors, '-o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'LineWidth', 2);
xlabel('c');
ylabel('Relative Error');
grid on;
box on;



function dydz = odes(z, y, c)
    U = y(1);
    V = y(2);

    a = 0.1;

    dUdz = V;
    dVdz = -c*V-U*(1-U)*(U-a);
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

    a = 0.1;

    du_dy = [0; (u(3:end) - u(1:end-2))/(2*dy); 0];
    d2u_dy2 = [0; (u(3:end) - 2*u(2:end-1) + u(1:end-2))/(dy^2); 0];
    
    f = alpha^(-1)*du_dy.*(p.c.*p.y)+alpha^(-2)*d2u_dy2+u.*(1-u).*(u-a); % RHS of the PDE
end

kappa_values = [0, 5, 20, 30, 50];

Y = 1;
N1 = 2e+5;
dy = Y/(N1-1);
y = linspace(0, Y, N1);

T = 100;
N2 = 1e+5;
dt = T/N2;

M = 1e+5;
L0 = 1;
L = L0;
alpha = M+L;

p.dy = dy;
p.dt = dt;
p.y = y(:);
p.L = L;
p.alpha = alpha;
p.M = M;

colors = ['r', 'g', 'b', 'y', 'm'];    % Colors for each kappa
lineStyles = {'-', '-', '-', '-', '-'};  % Line styles
%lineStyles = {'-', '--', '-.', ':', '-'};  % Line styles

figure; hold on;  % Prepare figure to hold all plots

for idx = 1:length(kappa_values)
    p.kappa = kappa_values(idx);  % Set current kappa value
    u = zeros(N1, 1);
    u(:, 1) = 1*(1-(y>1));  % Reset initial condition
    
    % RK4 time-stepping for current kappa value
    for n = 1:N2
        du_dy_at_y1 = (u(N1)-u(N1-1))/p.dy;
        p.c = -p.kappa/p.alpha*du_dy_at_y1;
        
        u = RK4total(p, u);
        
        p.L = p.L+p.c*p.dt;
        p.alpha = p.M+p.L;

        u(1) = u(2);
        u(end) = 0;
    end

    x = -p.M+p.y.*(p.L+p.M);
    U = u(:, end);
    z = x-p.c*T;
    I = find(U<1e-6, 1, 'first');
    z0 = z(I);

    plot(x-z0, U, 'LineWidth', 2, 'Color', colors(idx), 'LineStyle', lineStyles{idx}); % Plot with specific color and line style
    pause(0.01)
end


xlim([-20, 0]);
ylim([0, 1]);
xlabel('z');
ylabel('U');
legend(arrayfun(@(k) ['\kappa = ' num2str(k)], kappa_values, 'UniformOutput', false));  % Add legend
box on;
hold off;


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

    a = 0.1; % Example a value

    du_dy = [0; (u(3:end) - u(1:end-2))/(2*dy); 0];
    d2u_dy2 = [0; (u(3:end) - 2*u(2:end-1) + u(1:end-2))/(dy^2); 0];
    
    f = alpha^(-1)*du_dy.*(p.c.*p.y)+alpha^(-2)*d2u_dy2+u.*(1-u).*(u-a); % RHS of the PDE
end

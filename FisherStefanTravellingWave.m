% PDE
kappa = 10;

Y = 1; % Y domain corresponds to the transformed X domain
N1 = 100000;
dy = Y/(N1-1);
y = linspace(0, Y, N1);

T = 10;
N2 = 10000;
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

c = p.c;
disp(['Estimated wave speed c = ', num2str(c)]);


x = -p.M+p.y.*(p.L+p.M);

z = x-c*T;

[~, I_PDE_half] = min(abs(u(:, end)-0.5)); % Find the closest value to U = 1/2
z_PDE_half = z(I_PDE_half); % Find its corresponding z
disp(z_PDE_half)



% U,V system
%c is calculated above
z1 = -1e+13;
z2 = 0;

epsilon = 1e-8;
U0 = 1-epsilon;
V0 = -epsilon;
y0 = [U0; V0];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[Z, Y] = ode15s(@(z, y) odes(z, y, c), [z1, z2], y0, options);

U = Y(:,1);
I = find(U<1e-6, 1, 'first');
z0 = Z(I);

Z_part = Z(1:I);
U_part = U(1:I);



kappa_calculated = FisherStefanPhasePlane(c);
disp(['Calculated kappa based on PhasePlane: ', num2str(kappa_calculated)]);



[~, I_UV_half] = min(abs(U-0.5)); % Find the closest value to U = 1/2
z_UV_half = Z(I_UV_half)-z0; % Find its corresponding z
disp(z_UV_half)

z_offset = z_UV_half - z_PDE_half;

disp(z_offset)



plot(z+z_offset, u(:, end), '-r', 'LineWidth', 2);
hold on;

plot(Z_part-z0, U_part, '--b', 'LineWidth', 2); % Plot the stable region

Zextend = linspace(-40, -20, 100);
Uextend = ones(100, 1);
plot(Zextend, Uextend, '-r', 'LineWidth', 2);

xlabel('z');
ylabel('U');
xlim([-20, 0]);
ylim([0,1]);
hold off;



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

    du_dy = [0; (u(3:end) - u(1:end-2))/(2*dy); 0];
    d2u_dy2 = [0; (u(3:end) - 2*u(2:end-1) + u(1:end-2))/(dy^2); 0];
    
    f = alpha^(-1)*du_dy.*(p.c.*p.y)+alpha^(-2)*d2u_dy2+u.*(1-u); % RHS of the PDE
end

function kappa = AlleeStefanPhasePlane(c)
    a = 0.1;
    ode = @(U, V) -c-U*(1-U)*(U-a)/V;

    epsilon = 1e-8;
    u = [1, 0]; % We integrate from U=1 to U=0
    V1 = -epsilon; % When U=1, V=0
    
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [U, V] = ode15s(ode, u, V1, options); % Solve the ODE
    
    V0 = interp1(U, V(:,1), 0, 'linear', 'extrap'); % The value of V at U=0
    kappa = -c/V0;
    
    % disp(['Calculated kappa: ', num2str(kappa)]);
end
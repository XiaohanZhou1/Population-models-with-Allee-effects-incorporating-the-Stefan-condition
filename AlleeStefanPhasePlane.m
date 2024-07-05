% Input c, solve the ODE system to obtain kappa
function kappa = AEPhasePlane(c, a)
    z1 = 0;
    z2 = 1e+13;

    epsilon = 1e-8;
    U0 = 1-epsilon;
    V0 = -epsilon;
    y0 = [U0; V0];
    
    options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
    [Z, Y] = ode15s(@(z, y) odes(z, y, c, a), [z1, z2], y0, options);

    % Extract U and V from the solution Y
    U = Y(:,1);
    V = Y(:,2);
    
    % Find the index I where U is closest to zero
    [~, I] = min(abs(U));
    
    % Check if the corresponding V is close enough to be considered as V at U=0
    if abs(U(I)) < 1e-3 %1e-1
        % Calculate kappa using the formula kappa = -c/V0
        V0_at_U0 = V(I);  % V value when U is closest to zero
        kappa = -c / V0_at_U0;
        disp(['Calculated kappa = ', num2str(kappa)]);
    else
        disp('U=0 was not found in the calculated solution.');
        kappa = NaN;
    end
end

function dydz = odes(z, y, c, a)
    U = y(1);
    V = y(2);

    dUdz = V;
    dVdz = -c*V - U*(1-U)*(U-a);
    dydz = [dUdz; dVdz];
end

%cValues = [0.28, 0.50, 0.80];
%cValues = [0.56, 0.60, 0.80];
cValues = [0.54, 0.595];
z1 = -1e+13;
z2 = 0;
epsilon = 1e-8;
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
figure;
hold on;

% Store the points to be marked
markedPoints = [];
curvePoints = [];

for c = cValues
    y0 = [1-epsilon; -epsilon];

    [Z, Y] = ode15s(@(z, y) odes(z, y, c), [z1, z2], y0, options);

    U = Y(:,1);
    V = Y(:,2);

    plot(U, V, 'DisplayName', ['c = ', num2str(c)], 'LineWidth', 2);

    % Find the intersection of the curve with the y-axis (V-axis)
    foundIntersection = false;
    for i = 1:length(U)-1
        if U(i) * U(i+1) <= 0  % Check if the sign changes indicating an intersection
            % Linear interpolation to find the intersection point
            intersectV = V(i) - U(i) * (V(i+1) - V(i)) / (U(i+1) - U(i));
            curvePoints = [curvePoints; 0, intersectV];
            foundIntersection = true;
            break;
        end
    end

    if ~foundIntersection
        disp(['No intersection found for c = ', num2str(c)]);
    end
end

% Add the point (1, 0)
markedPoints = [markedPoints; 1, 0];

% Mark the specified points with circles and add legends with coordinates
plot(markedPoints(:,1), markedPoints(:,2), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', ...
    'DisplayName', '(1, 0)');

% Draw the x and y axes
plot([-0.2, 1.2], [0, 0], 'k', 'LineWidth', 0.5);
plot([0, 0], [-0.2, 0.1], 'k', 'LineWidth', 0.5);

xlabel('U');
ylabel('V');
box on;

if ~isempty(curvePoints)
    % Update the legend for the second point
    secondPointLegend = ['(0, ', num2str(curvePoints(1, 2), '%.4f'), ')'];
    plot(curvePoints(:,1), curvePoints(:,2), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', ...
        'DisplayName', secondPointLegend);
    legend('show');
else
    legend('show');
    disp('No curve points were found across all values of c.');
end

% Remove grid lines but keep x and y axes
ax = gca;
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
ax.XColor = 'k';
ax.YColor = 'k';
ax.XGrid = 'off';
ax.YGrid = 'off';

% Set xlim and ylim for better visualization
xlim([-0.2 1.2]);
ylim([-0.2 0.1]);

hold off;
  
function dydz = odes(z, y, c)
    U = y(1);
    V = y(2);

    a = 0.08;

    dUdz = V;
    dVdz = -c*V - U*(1-U)*(U-a);
    dydz = [dUdz; dVdz];
end








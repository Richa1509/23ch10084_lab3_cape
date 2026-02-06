%% Part (c): Solving with MATLAB bvp4c
clear; clc;

% --- Given Parameters ---
L = 10;
T0 = 300;
TL = 400;
Tinf = 200;
alpha = 0.05;
beta = 2.7e-9;

% 1. Create an initial guess (mesh and temperature)
% bvpinit(mesh, initial_guess_vector)
% We guess a linear temperature profile and a constant slope
solinit = bvpinit(linspace(0, L, 10), [350; 10]);

% 2. Solve the BVP
% @fin_ode: The system of first-order ODEs
% @fin_bc:  The boundary condition residuals
sol = bvp4c(@(x,y) fin_ode(x, y, alpha, beta, Tinf), ...
            @(ya,yb) fin_bc(ya, yb, T0, TL), ...
            solinit);

% 3. Evaluate the solution and Print Results
x_eval = linspace(0, L, 11);
T_sol = deval(sol, x_eval);

fprintf('--- Part (c): bvp4c Results ---\n');
fprintf('%-10s | %-10s\n', 'x', 'T(x)');
fprintf('---------------------------\n');
for i = 1:length(x_eval)
    fprintf('%-10.1f | %-10.4f\n', x_eval(i), T_sol(1,i));
end

% 4. Plotting
plot(sol.x, sol.y(1,:), 'b-o', 'LineWidth', 1.5);
xlabel('Length of fin (x)');
ylabel('Temperature T(x)');
title('Temperature Distribution along Fin (bvp4c)');
grid on;

%% --- ODE Function ---
function dydx = fin_ode(x, y, a, b, Tinf)
    % y(1) = T
    % y(2) = dT/dx
    % Convert 2nd order to 1st order:
    % dT/dx = y(2)
    % d2T/dx2 = alpha*(T - Tinf) + beta*(T^4 - Tinf^4)
    dydx = [y(2);
            a*(y(1) - Tinf) + b*(y(1)^4 - Tinf^4)];
end

%% --- Boundary Condition Function ---
function res = fin_bc(ya, yb, T0, TL)
    % ya: state at x=0, yb: state at x=L
    % Residuals should be zero when conditions are met
    res = [ya(1) - T0;   % T(0) = T0
           yb(1) - TL];  % T(L) = TL
end
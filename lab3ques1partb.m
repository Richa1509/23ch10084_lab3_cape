%% Problem 1B: Stiff Differential Equations Comparison
clear; clc; close all;
% Initial conditions and time span
y0 = [52.29; 83.82]; % [cite: 34]
tspan = [0 10]; % Define a reasonable time span to observe the decay
% --- 1. Solve using ode45 (Explicit Runge-Kutta) ---
tic; % Start timer
[t45, y45] = ode45(@stiff_system, tspan, y0);
time45 = toc; % End timer
% --- 2. Solve using ode15s (Stiff Solver) ---
tic; % Start timer
[t15, y15] = ode15s(@stiff_system, tspan, y0);
time15 = toc; % End timer
% --- 3. Visualization and Comparison ---
figure;
subplot(2,1,1);
plot(t45, y45(:,1), 'r', t15, y15(:,1), 'b--');
title('Comparison of y_1(t)'); ylabel('y_1'); legend('ode45', 'ode15s');
subplot(2,1,2);
plot(t45, y45(:,2), 'r', t15, y15(:,2), 'b--');
title('Comparison of y_2(t)'); ylabel('y_2'); xlabel('Time (t)'); legend('ode45', 'ode15s');
% --- 4. Print Efficiency Report ---
fprintf('--- Efficiency Report ---\n');
fprintf('ode45 Steps: %d | Time: %.4f seconds\n', length(t45), time45);
fprintf('ode15s Steps: %d | Time: %.4f seconds\n', length(t15), time15);
% --- Derivative Function ---
function dydt = stiff_system(~, y)
    y1 = y(1);
    y2 = y(2);
    
    dy1dt = -5*y1 + 3*y2;       % 
    dy2dt = 100*y1 - 301*y2;    % 
    
    dydt = [dy1dt; dy2dt];
end
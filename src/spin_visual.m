clc;
clear;

% Load constants and data
[h, hbar, kB, gyro, gyro_bar, I, Ns, mI, B0, Ts, w0, t] = constants();
load('spin_sys.mat', 'mu_x_all', 'mu_y_all', 'mu_z_all', 'mu_xt_all', 'mu_yt_all', 'mu_zt_all');
load('bulk_magnetization.mat', 'M_magn');

%% Plot 1: Individual Magnetic Moments and Their Precession

% Visualization setup for magnetic moments and their precession
figure;
hold on;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Magnetic Moments and Precession');
axis equal;
view(3);

% Colors for different vectors
colors = ['r', 'b'];  % Red for +1/2, Blue for -1/2
scaling_factor = 1e17;

% Plot the individual magnetic moment vectors and their precession
for i = 1:length(mI)
    % Apply scaling to the loaded data
    mu_x_scaled = mu_x_all(i, :) * scaling_factor;
    mu_y_scaled = mu_y_all(i, :) * scaling_factor;
    mu_z_scaled = mu_z_all(i, :) * scaling_factor;
    mu_xt_scaled = mu_xt_all(i, :) * scaling_factor;
    mu_yt_scaled = mu_yt_all(i, :) * scaling_factor;
    mu_zt_scaled = mu_zt_all(i, :) * scaling_factor;

    % Plot the non-scaled vector at the initial position
    quiver3(0, 0, 0, mu_x_scaled, mu_y_scaled, mu_z_scaled, colors(i), 'LineWidth', 2, 'MaxHeadSize', 0.5);
    plot3(mu_x_scaled, mu_y_scaled, mu_z_scaled, 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(i));

    % Plot the precession trajectory over time
    plot3(mu_xt_scaled, mu_yt_scaled, mu_zt_scaled, colors(i), 'LineWidth', 2);
end

% Add a legend for clarity
legend('mI = +1/2 (Initial)', 'mI = -1/2 (Initial)', 'mI = +1/2 (Precession)', 'mI = -1/2 (Precession)');

hold off;

%% Plot 2: Bulk Magnetization Vector

% Visualization of the bulk magnetization vector
figure;
quiver3(0, 0, 0, 0, 0, M_magn, 'k', 'LineWidth', 3, 'MaxHeadSize', 0.5);
hold on;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Bulk Magnetization Vector');
axis equal;
view(3);

% Add a dot at the end of the vector to indicate the endpoint
plot3(0, 0, M_magn, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');

legend('Bulk Magnetization');
hold off;

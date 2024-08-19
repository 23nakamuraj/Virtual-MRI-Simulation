
% Initialization
t = linspace(0, 1, 100); % time (s)

T1 = 0.8;  T2 = 0.1; % Relaxation time (s)
M0 = 1; % Initial magnetization

M_start = [M0, 0, 0].'; % Start along Mx (assume initial magnetization along x-axis)

Mall = zeros(3, length(t)); % Pre-allocate matrix to store magnetization vectors

% Simulate Bloch relaxation over time
for It = 1:length(t)
    Mall(:, It) = bloch_relax(M_start, t(It), M0, T1, T2);
end

% Plot the results
figure;
plot(t, Mall(1,:), 'r', 'DisplayName', 'Mx');
hold on;
plot(t, Mall(2,:), 'g', 'DisplayName', 'My');
plot(t, Mall(3,:), 'b', 'DisplayName', 'Mz');
xlabel('Time (s)');
ylabel('Magnetization');
legend;
title('Bloch Relaxation');
grid on;



function M = bloch_relax(M_start, t, M0, T1, T2)
    % T1 relaxation (longitudinal relaxation)
    Mz = M0 * (1 - exp(-t/T1));
    
    % T2 relaxation (transverse relaxation)
    Mxy = M_start(1:2) * exp(-t/T2);
    
    % Combining longitudinal and transverse components
    M = [Mxy; Mz];
end
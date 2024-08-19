clc;
clear;

[h, hbar, kB, gyro, gyro_bar, I, Ns, mI, B0, Ts, w0, t] = constants();

function [E_up, E_down] = spin_up_down(mI, gyro, hbar,B0)
    E_up = -mI(1) * gyro * hbar * B0;
    E_down = -mI(2) * gyro * hbar * B0;
end

function E_delta = spin_energy_difference(E_up, E_down)
    E_delta = E_down - E_up;
end

function [N_ratio, N_diff, N_up, N_down] = spin_population(E_delta,gyro,hbar,kB,Ts,Ns,B0)
    N_ratio = exp(-E_delta/(kB*Ts));
    N_diff = (Ns * gyro * hbar * B0)/(2 * kB * Ts);
    N_up = N_diff/(1-N_ratio);
    N_down = N_up / N_ratio;
end

function M = bulk_magnetization_vector(N_up, N_down, gyro, hbar)
    M = 0.5 * (N_up - N_down) * gyro * hbar;
end

function M_magn = bulk_magnetization_magnitude(gyro, hbar, B0, Ns, I, kB, Ts)
    M_magn = (gyro^2 * hbar^2 * B0 * Ns* I * (I + 1))/(3* kB * Ts); 
end
  
[E_up, E_down] = spin_up_down(mI, gyro, hbar, B0);
E_delta = spin_energy_difference(E_up, E_down);
[N_ratio, N_diff, N_up, N_down] = spin_population(E_delta,gyro,hbar,kB,Ts,Ns,B0);
M = bulk_magnetization_vector(N_up, N_down, gyro, hbar);
M_magn = bulk_magnetization_magnitude(gyro, hbar, B0, Ns, I, kB, Ts);

disp(['E_up is: ', num2str(E_up), ', E_down is: ', num2str(E_down),...
    ', E_delta is: ', num2str(E_delta)])
disp(['N_ratio = ', num2str(N_ratio), ', N_diff = ', num2str(N_up), ', N_up = ',...
    num2str(N_up), ', N_down = ', num2str(N_down)])
disp(['bulk magnetization vector = ', num2str(M)])
disp(['bulk magnetization magnitude= ', num2str(M_magn)])

save('bulk_magnetization.mat', 'M_magn');


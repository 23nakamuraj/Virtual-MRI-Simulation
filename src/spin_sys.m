clc;
clear;

% Load constants
[h, hbar, kB, gyro, gyro_bar, I, Ns, mI, B0, Ts, w0, t] = constants();

function mu = nuclear_magnetic_moment(gyro, hbar, I)
    mu = gyro * hbar * sqrt(I * (I + 1));
end

function mu_z = longitudinal_magnetic_moment(gyro, hbar, mI)
    mu_z = gyro * hbar * mI;
end

function mu_xy = transverse_magnetic_moment(gyro, hbar)
    mu_xy = gyro * hbar / sqrt(2);
end

function [mu_x, mu_y] = transverse_component(mu_xy)
    ang = rand * 2 * pi;
    mu_x = mu_xy * cos(ang);
    mu_y = mu_xy * sin(ang);
end

function theta = mu_orientation(mu, mu_z)
    theta = rad2deg(acos(mu_z / mu));
end

function [mu_xt, mu_yt, mu_zt] = nuclear_precession(mu_x, mu_y, mu_z, w0, t)
    mu_xt = mu_x * cos(w0*t) + mu_y * sin(w0*t);
    mu_yt = mu_y * cos(w0*t) - mu_x * sin(w0*t);
    mu_zt = mu_z * ones(size(t));  % Keep the z-component constant
end

% Initialize arrays to store results for both mI values
mu_x_all = [];
mu_y_all = [];
mu_z_all = [];
mu_xt_all = [];
mu_yt_all = [];
mu_zt_all = [];

% Calculate and store the spin system variables for each mI value
for i = 1:length(mI)
    mu = nuclear_magnetic_moment(gyro, hbar, I);
    mu_z = longitudinal_magnetic_moment(gyro, hbar, mI(i));
    mu_xy = transverse_magnetic_moment(gyro, hbar);
    [mu_x, mu_y] = transverse_component(mu_xy);
    [mu_xt, mu_yt, mu_zt] = nuclear_precession(mu_x, mu_y, mu_z, w0, t);

    % Store the results in the arrays
    mu_x_all = [mu_x_all; mu_x];
    mu_y_all = [mu_y_all; mu_y];
    mu_z_all = [mu_z_all; mu_z];
    mu_xt_all = [mu_xt_all; mu_xt];
    mu_yt_all = [mu_yt_all; mu_yt];
    mu_zt_all = [mu_zt_all; mu_zt];
end

% Save the results to a .mat file for later visualization
save('spin_sys.mat', 'mu_x_all', 'mu_y_all', 'mu_z_all', 'mu_xt_all', 'mu_yt_all', 'mu_zt_all');

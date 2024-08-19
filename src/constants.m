function [h, hbar, kB, gyro, gyro_bar, I, Ns, mI, B0, Ts, w0,t] = constants()
    h = 6.626e-34;        % Planck's constant in J*s
    hbar = h / (2 * pi);  % Reduced Planck's constant
    kB = 1.38e-23;       % Boltzmann constant in J/K
    gyro = 2.675e8;            % Gyromagnetic ratio for protons in rad/s/T
    gyro_bar = gyro / (2*pi); % Gyromagnetic ratio in Hz/T
    I = 0.5;  % Spin quantum number for protons
    Ns = 1e23;       % Number of spins
    mI = [+1/2, -1/2];  % Magnetic quantum number values for a spin-1/2 system (hydrogen)
    B0 = 1.5;        % External magnetic field in Tesla
    Ts = 310;        % Temperature in Kelvin
    w0 = gyro * B0; % Larmor frequency
    t = linspace(0, 1e-3, 1000);  % Time vector for one period
end


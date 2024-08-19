# Constants
h = 6.626e-34  # Planck's constant in J*s
hbar = h / (2 * np.pi)  # Reduced Planck's constant
k_B = 1.38e-23  # Boltzmann constant in J/K

# Variable constants
gamma_proton = 2.675e8
gamma_bar = gamma_proton/(2*np.pi)

# Gyromagnetic ratio rad/s/T
I_proton = 0.5  # Spin quantum number for protons
Ns = 1e23  # Number of spins
m_I_values = [+1/2, -1/2] # Magnetic quantum number values for a spin-1/2 system (hydrogen)

# MRI system variables
B0 = 1.5  # External magnetic field in Tesla
Ts = 310  # Temperature in Kelvin

# Larmor frequency
larmor = gamma_proton * B0

# Time array for simulation
t = np.linspace(0, 1e-3, 1000)  # 1 millisecond simulation

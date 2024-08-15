import numpy as np
import matplotlib.pyplot as plt

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




def nuclear_magnetic_moment(gamma, I):
    """
    Calculate the magnitude of the nuclear magnetic moment.
    
    Parameters:
    gamma (float): Gyromagnetic ratio in rad/s/T.
    I (float): Nuclear spin quantum number.
    
    Returns:
    float: Magnitude of the nuclear magnetic moment in Joules/Tesla.
    """
    mu = gamma * hbar * np.sqrt(I * (I + 1))
    return mu

def longitudinal_magnetic_moment(gamma, m_I):
    """
    Calculate the z-component of the magnetic moment.
    
    Parameters:
    gamma (float): Gyromagnetic ratio in rad/s/T.
    I (float): Nuclear spin quantum number.
    m_I (float): Magnetic quantum number.
    
    Returns:
    float: Z-component of the magnetic moment in Joules/Tesla.
    """
    mu_z = gamma * hbar * m_I
    return mu_z

def mu_orientation(mu_z, mu):
    """
    Calculate the angle theta between mu and B0. 
    We know that mu has two orientations by m1, the magnetic quantum number which takes (2I + 1) possible orientations, where I = 1/2.

    Parameters: 
    mu_z (float): Z-component of the magnetic moment in Joules/Tesla
    mu (float): Magnitude of the nuclear magnetic moment in Joules/Tesla

    Returns:
    float: Angle between mu_z and mu in degrees.
    """
    theta = np.degrees(np.arccos(mu_z/mu))
    return theta

def transverse_magnetic_moment(gamma):
    """
    Calculate the magnitude of the transverse magnetic moment (mu_x, mu_y).
    
    Parameters:
    gamma (float): Gyromagnetic ratio in rad/s/T.
    
    Returns:
    mu_xy (float): Magnitude of the transverse magnetic moment in Joules per Tesla.
    """
    mu_xy = gamma * hbar / np.sqrt(2)
    return mu_xy

def mu_xy_components(mu_xy):
    """
    Calculate the individual transverse components of mu (mu_x and mu_y) given that the direction of transverse component mu_xy remains random.
    In this case the angle is a random variable uniformly distributed over [0,2pi).
    Parameters: 
    mu_xy (float): Magnitude of the transverse magnetic moment in Joules per Tesla.

    Returns: 
    mu_x (float): Orientation of mu in the x-direction
    mu_y (float): Orientation of mu in the y-direction
    """
    angle = np.random.uniform(0, 2 * np.pi)
    mu_x = mu_xy * np.cos(angle)
    mu_y = mu_xy * np.sin(angle)
    return mu_x, mu_y


def nuclear_precession(mu_x_0, mu_y_0, mu_z_0):
    """
    Calculates the equation of motion for isolated spins in the classical treatment.
    LaTeX mode: $$\frac{\vec{\mu}}{dt} = \gamma\vec{mu}\times B_0\vec{k}$$

    Parameters: 
    mu_x_0 (float): Orientation of mu in the x-direction
    mu_y_0 (float): Orientation of mu in the y-direction
    mu_z_0 (float): Orientation of mu in the z-direction

    Results:
    mu_x_t (float): Time dependent orientation of mu in the x-direction
    mu_y_t (float): Time dependent orientation of mu in the y-direction
    mu_y_t (float): Time dependent orientation of mu in the z-direction
    """
    mu_x_t = mu_x_0 * np.cos(larmor * t) + mu_y_0 * np.sin(larmor * t)
    mu_y_t = mu_y_0 * np.cos(larmor * t) - mu_x_0 * np.sin(larmor * t)
    mu_z_t = mu_z_0 * np.ones_like(t)
    return mu_x_t, mu_y_t, mu_z_t

def rotation_matrix(mu_x_0, mu_y_0, mu_z_0):
    # Initial magnetic moment vector (assuming initially aligned along the x-axis)
    mu_initial = np.array([mu_x_0, mu_y_0, mu_z_0])
    
    # Precession using rotation matrix R_z(alpha) where alpha = larmor * t
    mu_x_t = []
    mu_y_t = []
    mu_z_t = []
    
    for time in t:
        alpha = larmor * time
        R_z = np.array([
            [np.cos(alpha), np.sin(alpha), 0],
            [-np.sin(alpha), np.cos(alpha), 0],
            [0, 0, 1]
        ])
        mu_t = np.dot(R_z, mu_initial)
        mu_x_t.append(mu_t[0])
        mu_y_t.append(mu_t[1])
        mu_z_t.append(mu_t[2])
    
    # Convert lists to numpy arrays for plotting
    mu_x_t_rot = np.array(mu_x_t)
    mu_y_t_rot = np.array(mu_y_t)
    mu_z_t_rot = np.array(mu_z_t)
    return mu_x_t_rot, mu_y_t_rot, mu_z_t_rot

mu_proton = nuclear_magnetic_moment(gamma_proton, I_proton)

# Calculate the z-component of the magnetic moment for both m_I values
mu_z_values = [longitudinal_magnetic_moment(gamma_proton, m_I) for m_I in m_I_values]

# Calculate the transverse magnetic moment
mu_xy = transverse_magnetic_moment(gamma_proton)

# Get the random components mu_x and mu_y
mu_x_0, mu_y_0 = mu_xy_components(mu_xy)
# Assume an initial z-component
mu_z_0 = longitudinal_magnetic_moment(gamma_proton, m_I_values[0])

# Calculate precession using second order differential equation
mu_x_t, mu_y_t, mu_z_t = nuclear_precession(mu_x_0, mu_y_0, mu_z_0)

# Calculate precession using rotational matrix
mu_x_t_rot, mu_y_t_rot, mu_z_t_rot = rotation_matrix(mu_x_0, mu_y_0, mu_z_0)

# Calculate the magnitude of the transverse magnetic moment from its components
magnitude_mu_xy = np.sqrt(mu_x_0**2 + mu_y_0**2)

# Calculate the orientation angles theta
theta_values = mu_orientation(mu_z_values, mu_proton)


# Display the results
print(f"mu_x: {mu_x_0:.3e} J/T, mu_y: {mu_y_0:.3e} J/T")
print(f'mu_z: {mu_z_values}')
print(f"Calculated|μ_xy|: {magnitude_mu_xy:.3e} J/T")
# Ensure it matches the initial calculated mu_xy magnitude
print(f"Expected |μ_xy|: {mu_xy:.3e} J/T")

# Display results
for m_I, mu_z, theta in zip(m_I_values, mu_z_values, theta_values):
    print(f"For m_I = {m_I}: μ_z = {mu_z:.3e} J/T, θ = {theta:.2f} degrees")
    

print(f"Nuclear magnetic moment (mu) for protons: {mu_proton:.3e} J/T")


# Plotting the precession
plt.figure(figsize=(10, 6))
plt.plot(t * 1e3, mu_x_t, label='μ_x(t)', color='blue')
plt.plot(t * 1e3, mu_y_t, label='μ_y(t)', color='red')
plt.plot(t * 1e3, mu_z_t, label='μ_z(t)', color='yellow')
plt.title('Nuclear Precession Over Time')
plt.xlabel('Time (ms)')
plt.ylabel('Magnetic Moment (Arbitrary Units)')
plt.legend()
plt.grid(True)
plt.show()

# Plotting the precession
plt.figure(figsize=(10, 6))
plt.plot(t * 1e3, mu_x_t_rot, label='μ_x(t)', color='blue')
plt.plot(t * 1e3, mu_y_t_rot, label='μ_y(t)', color='red')
plt.plot(t * 1e3, mu_z_t_rot, label='μ_z(t)', color='yellow')
plt.title('Nuclear Precession Using Rotation Matrix R_z(α)')
plt.xlabel('Time (ms)')
plt.ylabel('Magnetic Moment (Arbitrary Units)')
plt.legend()
plt.grid(True)
plt.show()






def spin_up_down(m_I):
    """
    Applying quantum theory E = -mu * B0 for pointing-up spins(m=1/2)

    Parameters:
    m_I_values (list): Contains the magnetic quantum numbers mI = 1/2, -1/2
    
    Result:
    E_up (float): Energy of interaction for point up spins (m=1/2)
    E_down (float): Energy of interaction for point down spins (m=-1/2)
    """
    E_up = -m_I_values[0] * gamma_proton * hbar * B0
    E_down = m_I_values[0] * gamma_proton * hbar * B0
    return E_up, E_down

def spin_energy_diff(E_up, E_down):
    """
    Calculates the energy difference between the two spin states.
    The following equation calculates the same: gamma_proton * hbar * B0

    Result:
    E_delta (float): Energy difference between E_up and E_down.
    """
    E_delta = E_down - E_up
    return E_delta

def spin_population(E_delta):
    """
    Calculate the population ratio and difference between spin-up and spin-down states.
    
    Parameters:
    E_delta (float): Energy difference between E_up and E_down.
    
    Returns:
    N_up (float): Number of spins in the up state.
    N_down (float): Number of spins in the down state.
    """
    # Population ratio (Boltzmann distribution)
    N_ratio = np.exp(-E_delta / (k_B * Ts))
    
    # Total population difference
    N_diff = Ns * gamma_proton * hbar * B0 / (2 * k_B * Ts)
    
    # Solve for N_up and N_down
    N_up = N_diff / (1 - N_ratio)
    N_down = N_up / N_ratio
    return N_up, N_down, N_ratio, N_diff

def bulk_magnetization(N_up, N_down, m_I_values):
    """
    Calculates the bulk magnetization vector from spin system using eq below.
    M = M_x{i} + M_y{j} + M_z{k}
    The bulk magnetization vector points exactly along the positive direction of the z-axis at equilibrium.

    Parameters:
    

    Returns:
    
    """
    M = 0.5*(N_up - N_down) * gamma_proton * hbar
    return M

def bulk_magnetization_magnitude(M):
    M_magn = (gamma_proton**2 * hbar**2 * B0 * Ns* I_proton * (I_proton + 1))/(3* k_B * Ts) 
    return M_magn

# Calculate the energy levels
E_up, E_down = spin_up_down(m_I_values)

# Calculate the energy difference
E_delta = spin_energy_diff(E_up, E_down)

# Calculate the spin populations
N_up, N_down, N_ratio, N_diff = spin_population(E_delta)

M  = bulk_magnetization(N_up, N_down, m_I_values)
M_magn = bulk_magnetization_magnitude(M)
print(f'Bulk magnetization is {M} and its magnitude is:{M_magn}')
print(f'E_up is: {E_up}\nE_down is: {E_down} \nDifference is:  {E_delta}')
print(f'The number of spins is {Ns} with spin up {N_up} and spin down {N_down}')
print(f'The ratio between spins is {N_ratio} and diff {N_diff}')

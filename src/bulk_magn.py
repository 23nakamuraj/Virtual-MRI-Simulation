from constant import *
import numpy as np
import matplotlib.pyplot as plt

def spin_up_down():
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

def bulk_magnetization(N_up, N_down):
    """
    Calculates the bulk magnetization vector from spin system using eq below.
    M = M_x{i} + M_y{j} + M_z{k}
    The bulk magnetization vector points exactly along the positive direction of the z-axis at equilibrium.

    Parameters:
    

    Returns:
    
    """
    M = 0.5*(N_up - N_down) * gamma_proton * hbar
    return M

def bulk_magnetization_magnitude():
    M_magn = (gamma_proton**2 * hbar**2 * B0 * Ns* I_proton * (I_proton + 1))/(3* k_B * Ts) 
    return M_magn

# Calculate the energy levels
E_up, E_down = spin_up_down()

# Calculate the energy difference
E_delta = spin_energy_diff(E_up, E_down)

# Calculate the spin populations
N_up, N_down, N_ratio, N_diff = spin_population(E_delta)

M  = bulk_magnetization(N_up, N_down)
M_magn = bulk_magnetization_magnitude()
print(f'Bulk magnetization is {M} and its magnitude is:{M_magn}')
print(f'E_up is: {E_up}\nE_down is: {E_down} \nDifference is:  {E_delta}')
print(f'The number of spins is {Ns} with spin up {N_up} and spin down {N_down}')
print(f'The ratio between spins is {N_ratio} and diff {N_diff}')


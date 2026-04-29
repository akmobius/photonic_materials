import numpy as np
import matplotlib.pyplot as plt

def meta_atom_tmm(wavelength, thickness, n_core, n_substrate):
    """
    Computes transmission and phase delay for a meta-atom using FP logic.
    Ref: Lecture 13, TMM Examples (pg 7-8)
    """
    k0 = 2 * np.pi / wavelength
    
    # 1. Define Admittances (Normal Incidence for simplicity)
    # Ref: Source [1234]
    eta0 = 376.73  # Vacuum impedance
    Y0 = 1.0 / eta0             # Air
    Ys = n_substrate / eta0     # Substrate
    Yc = n_core / eta0          # Meta-atom core
    
    # 2. Interface Matrices (D) - Entrance and Exit
    # Ref: Source [1283, 3687]
    D01 = 0.5 * np.array([[1 + Yc/Y0, 1 - Yc/Y0],
                          [1 - Yc/Y0, 1 + Yc/Y0]])
    
    D1s = 0.5 * np.array([[1 + Ys/Yc, 1 - Ys/Yc],
                          [1 - Ys/Yc, 1 + Ys/Yc]])
    
    # 3. Propagation Matrix (P) through height h
    # Ref: Source [1303, 3689]
    phi = n_core * k0 * thickness
    P1 = np.array([[np.exp(1j * phi), 0],
                   [0, np.exp(-1j * phi)]])
    
    # 4. System Transfer Matrix (M)
    # Ref: Source [1320, 3693]
    M = D01 @ P1 @ D1s
    
    # 5. Extract Complex Transmission t = 1/M11
    # Ref: Source [1325, 3694]
    t = 1.0 / M[0, 0]
    
    transmission = np.abs(t)**2
    phase = np.angle(t) # in radians
    
    return transmission, phase

# --- Simulation Parameters ---
wavelength = 0.785 # microns (HeNe Laser) [cite: 86]
n_Si = 3.48        # Silicon meta-atom [cite: 3190]
n_glass = 1.45     # Silica substrate [cite: 3190, 1654]
heights = np.linspace(0.1, 1.0, 500) # microns

# --- Run Sweep ---
T_list = []
P_list = []

for h in heights:
    trans, ph = meta_atom_tmm(wavelength, h, n_Si, n_glass)
    T_list.append(trans)
    P_list.append(ph)

# --- Plotting ---
fig, ax1 = plt.subplots(figsize=(10, 6))

color = 'tab:blue'
ax1.set_xlabel('Meta-atom Height (microns)')
ax1.set_ylabel('Transmission Intensity', color=color)
ax1.plot(heights, T_list, color=color, linewidth=2, label='Transmission')
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid(True, linestyle='--', alpha=0.6)

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel('Phase Delay (radians)', color=color)
# Unwrapping phase to see continuous trend
ax2.plot(heights, np.unwrap(P_list), color=color, linestyle='--', label='Phase')
ax2.tick_params(axis='y', labelcolor=color)

plt.title(f'Fabry-Pérot Meta-atom Characterization (λ={wavelength} microns)')
fig.tight_layout()
plt.show()
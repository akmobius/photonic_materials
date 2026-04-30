import numpy as np
import matplotlib.pyplot as plt

def get_meta_atom_results(wavelength, h, n_core, n_sub, pol='TE'):
    """
    Computes T and Phase using your original TMM logic structure.
    Adjusted for Polarization and 785nm NIR parameters.
    """
    
    """Computes T and Phase using TMM logic"""
    k0 = 2 * np.pi / wavelength
    eta0 = 376.73
    Y0, Yc, Ys = 1.0/eta0, n_core/eta0, n_sub/eta0
    
    # Matching matrices D 
    D01 = 0.5 * np.array([[1 + Yc/Y0, 1 - Yc/Y0], [1 - Yc/Y0, 1 + Yc/Y0]])
    D1s = 0.5 * np.array([[1 + Ys/Yc, 1 - Ys/Yc], [1 - Ys/Yc, 1 + Ys/Yc]])
    
    # Propagation matrix P through height h
    phi = n_core * k0 * h
    P1 = np.array([[np.exp(1j * phi), 0], [0, np.exp(-1j * phi)]])
    
    # System Matrix M
    M = D01 @ P1 @ D1s
    
    # Transmission amplitude t = 1/M11 
    t = 1.0 / M[0, 0]
    return np.abs(t)**2, np.angle(t)

# --- Configuration for 785nm Laser ---
wl_design = 0.785  # microns
n_Si_785 = 3.71    # Silicon at 785nm
n_glass = 1.45     # Holder substrate
h_nominal = 0.8    # microns
h_error = h_nominal * 1.10  # 10% fabrication error study

# Range of effective indices representing varying meta-atom radii
n_eff_sweep = np.linspace(1.2, n_Si_785 - 0.5, 500)

# --- Calculations ---
T_nom, P_nom = [], []
T_err, P_err = [], []
T_tm, P_tm = [] , [] # For Anisotropy Study

for neff in n_eff_sweep:
    # 1. Nominal Study (TE/Standard)
    tn, pn = get_meta_atom_results(wl_design, h_nominal, neff, n_glass, 'TE')
    T_nom.append(tn); P_nom.append(pn)
    
    # 2. 10% Height Sensitivity Study
    te, pe = get_meta_atom_results(wl_design, h_error, neff, n_glass, 'TE')
    T_err.append(te); P_err.append(pe)
    
    # 3. Anisotropy Study (TM/p-pol)
    # Assumes meta-atom is non-cylindrical, providing a shifted effective index
    ttm, ptm = get_meta_atom_results(wl_design, h_nominal, neff * 0.95, n_glass, 'TM')
    T_tm.append(ttm); P_tm.append(ptm)

# --- Holder Spectral Dispersion (Kramers-Kronig logic) ---
# Modeling index shift over a broadband 700nm-900nm range
broadband_wl = np.linspace(0.700, 0.900, 200)
n_holder = 1.45 + (0.01 / broadband_wl**2) # Simplified dispersion formula

# --- Plotting ---
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Transmission Comparison (Sensitivity & Anisotropy)
axs[0, 0].plot(n_eff_sweep, T_nom, 'b', label='Nominal (h=0.8)')
axs[0, 0].plot(n_eff_sweep, T_err, 'r--', label='10% Error (h=0.88)')
axs[0, 0].set_title("Transmission Resonance Shifts")
axs[0, 0].set_xlabel("Effective Index (Radius Proxy)")
axs[0, 0].set_ylabel("Intensity")
axs[0, 0].legend()
axs[0, 0].grid(True, alpha=0.3)

# Plot 2: Phase Comparison
axs[0, 1].plot(n_eff_sweep, np.unwrap(P_nom), 'b', label='Nominal Phase')
axs[0, 1].plot(n_eff_sweep, np.unwrap(P_err), 'r--', label='Error Phase')
axs[0, 1].set_title("Phase Sensitivity Study")
axs[0, 1].set_xlabel("Effective Index")
axs[0, 1].set_ylabel("Phase (rad)")
axs[0, 1].legend()

# Plot 3: Anisotropy Study (Polarization Phase Difference)
axs[1, 0].plot(n_eff_sweep, np.unwrap(P_nom), 'g', label='TE-Phase')
axs[1, 0].plot(n_eff_sweep, np.unwrap(P_tm), 'm', label='TM-Phase')
axs[1, 0].set_title("Anisotropy: Polarization Splitting")
axs[1, 0].set_xlabel("Effective Index")
axs[1, 0].legend()

# Plot 4: Holder Dispersion (Imaging Impact)
axs[1, 1].plot(broadband_wl, n_holder, 'k')
axs[1, 1].set_title("Holder Dispersion (Broadband Performance)")
axs[1, 1].set_xlabel("Wavelength (microns)")
axs[1, 1].set_ylabel("Substrate Index (n)")

plt.tight_layout()
plt.show()
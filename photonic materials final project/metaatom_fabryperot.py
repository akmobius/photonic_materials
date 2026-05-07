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

# --- 4. Dispersion Impact on Characterization Sweep ---
P_700, P_900 = [], []
T_700, T_900 = [], []
wl_700 = 0.700
n_sub_700 = 1.45 + (0.01 / wl_700**2)
wl_900 = 0.900
n_sub_900 = 1.45 + (0.01 / wl_900**2)

for neff in n_eff_sweep:
    t700, p700 = get_meta_atom_results(wl_700, h_nominal, neff, n_sub_700, 'TE')
    t900, p900 = get_meta_atom_results(wl_900, h_nominal, neff, n_sub_900, 'TE')
    P_700.append(p700)
    P_900.append(p900)
    T_700.append(t700)
    T_900.append(t900)

# --- Plotting ---
# Original Figure
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
axs[1, 1].set_title("Holder Dispersion Profile (Kramers-Kronig)")
axs[1, 1].set_xlabel("Wavelength (microns)")
axs[1, 1].set_ylabel("Substrate Index (n)")

plt.tight_layout()

# --- New Separate Figure for Dispersion Impact ---
fig2, (ax_t, ax_p) = plt.subplots(1, 2, figsize=(12, 5))

# Transmission with Dispersion
ax_t.plot(n_eff_sweep, T_700, 'c', label='700nm')
ax_t.plot(n_eff_sweep, T_nom, 'b', label='785nm')
ax_t.plot(n_eff_sweep, T_900, 'm', label='900nm')
ax_t.set_title("Dispersion Impact on Transmission")
ax_t.set_xlabel("Effective Index")
ax_t.set_ylabel("Transmission Intensity")
ax_t.legend()
ax_t.grid(True, alpha=0.3)

# Phase with Dispersion
ax_p.plot(n_eff_sweep, np.unwrap(P_700), 'c', label='700nm')
ax_p.plot(n_eff_sweep, np.unwrap(P_nom), 'b', label='785nm')
ax_p.plot(n_eff_sweep, np.unwrap(P_900), 'm', label='900nm')
ax_p.set_title("Dispersion Impact on Phase Library")
ax_p.set_xlabel("Effective Index")
ax_p.set_ylabel("Phase (rad)")
ax_p.legend()
ax_p.grid(True, alpha=0.3)

plt.tight_layout()

# --- Combined Phase and Transmission Sweep (Radius Proxy) ---
fig3 = plt.figure(figsize=(14, 5))

# Plot 1: Geometry (Top-Down Radial Profile)
ax_geo = fig3.add_subplot(1, 2, 1)
x = np.linspace(-0.5, 0.5, 1000)
sample_radii = [0.15, 0.25, 0.35]
colors = ['lightblue', 'skyblue', 'steelblue']
for r, c in zip(sample_radii, colors):
    n_profile_x = np.piecewise(x, [np.abs(x) <= r, np.abs(x) > r], [n_Si_785, 1.0])
    ax_geo.plot(x, n_profile_x, color=c, lw=2, label=f'Radius = {r} $\mu$m')
    ax_geo.fill_between(x, 1.0, n_profile_x, color=c, alpha=0.2)

ax_geo.set_title("Assumed Geometry (Top-Down Radial Profile)")
ax_geo.set_xlabel("Radial coordinate x (microns)")
ax_geo.set_ylabel("Refractive Index (n)")
ax_geo.set_ylim(0.8, n_Si_785 + 0.5)
ax_geo.annotate('Air (n=1)', xy=(-0.45, 1.1), fontweight='bold')
ax_geo.annotate(f'Si Pillar (n={n_Si_785})', xy=(-0.15, n_Si_785 + 0.1), fontweight='bold')
ax_geo.legend(loc='lower right')
ax_geo.grid(True, alpha=0.3)

# Plot 2: Combined Sweep
ax_res = fig3.add_subplot(1, 2, 2)
color1 = 'tab:blue'
ax_res.plot(n_eff_sweep, T_nom, color=color1, lw=2, label='Transmission')
ax_res.set_xlabel('Effective Index ($n_{eff}$ - proxy for Pillar Radius)')
ax_res.set_ylabel('Transmission Intensity', color=color1)
ax_res.tick_params(axis='y', labelcolor=color1)

ax_ph = ax_res.twinx()
color2 = 'tab:red'
ax_ph.plot(n_eff_sweep, np.unwrap(P_nom), color=color2, ls='--', label='Phase')
ax_ph.set_ylabel('Phase Delay (radians)', color=color2)
ax_ph.tick_params(axis='y', labelcolor=color2)
ax_res.set_title("Combined Characterization Sweep")

fig3.tight_layout()

plt.show()
import numpy as np
import matplotlib.pyplot as plt

def get_meta_atom_results(wavelength, h_range, n_core, n_sub):
    """Computes T and Phase using TMM logic [Source: 3692-3694]"""
    k0 = 2 * np.pi / wavelength
    eta0 = 376.73
    Y0, Yc, Ys = 1.0/eta0, n_core/eta0, n_sub/eta0
    
    # Matching matrices D [Source: 1283, 3687]
    D01 = 0.5 * np.array([[1 + Yc/Y0, 1 - Yc/Y0], [1 - Yc/Y0, 1 + Yc/Y0]])
    D1s = 0.5 * np.array([[1 + Ys/Yc, 1 - Ys/Yc], [1 - Ys/Yc, 1 + Ys/Yc]])
    
    trans_list, phase_list = [], []
    for h in h_range:
        # Propagation matrix P [Source: 1303, 3689]
        phi = n_core * k0 * h
        P1 = np.array([[np.exp(1j * phi), 0], [0, np.exp(-1j * phi)]])
        M = D01 @ P1 @ D1s
        t = 1.0 / M[0, 0] # [Source: 3694]
        trans_list.append(np.abs(t)**2)
        phase_list.append(np.angle(t))
    return np.array(trans_list), np.array(phase_list)

# --- Configuration ---
wl = 0.785
n_Si, n_glass = 3.48, 1.45
h_sweep = np.linspace(0.1, 1.0, 500)
sample_h = 0.6 # The specific height for geometry plot

# --- Calculations ---
T, P = get_meta_atom_results(wl, h_sweep, n_Si, n_glass)

# --- Plotting ---
fig = plt.figure(figsize=(12, 5))

# Plot 1: Geometry (Refractive Index Profile)
ax_geo = fig.add_subplot(1, 2, 1)
z = np.linspace(-0.2, sample_h + 0.2, 1000)
n_profile = np.piecewise(z, [z < 0, (z >= 0) & (z <= sample_h), z > sample_h], 
                         [1.0, n_Si, n_glass])

ax_geo.plot(z, n_profile, 'black', lw=2)
ax_geo.fill_between(z, 0, n_profile, color='orange', alpha=0.2)
ax_geo.set_title(f"Assumed Geometry (h={sample_h} microns)")
ax_geo.set_xlabel("z-coordinate (microns)")
ax_geo.set_ylabel("Refractive Index (n)")
ax_geo.annotate('Air (n=1)', xy=(-0.1, 1.1), fontweight='bold')
ax_geo.annotate(f'Si Pillar (n={n_Si})', xy=(sample_h/4, n_Si-0.3), fontweight='bold')
ax_geo.annotate(f'Substrate (n={n_glass})', xy=(sample_h+0.02, n_glass+0.1), fontweight='bold')
ax_geo.grid(True, alpha=0.3)

# Plot 2: Phase and Transmission Sweep
ax_res = fig.add_subplot(1, 2, 2)
color1 = 'tab:blue'
ax_res.plot(h_sweep, T, color=color1, lw=2, label='Transmission')
ax_res.set_xlabel('Pillar Height (microns)')
ax_res.set_ylabel('Transmission Intensity', color=color1)
ax_res.tick_params(axis='y', labelcolor=color1)

ax_ph = ax_res.twinx()
color2 = 'tab:red'
ax_ph.plot(h_sweep, np.unwrap(P), color=color2, ls='--', label='Phase')
ax_ph.set_ylabel('Phase Delay (radians)', color=color2)
ax_ph.tick_params(axis='y', labelcolor=color2)

plt.title("FP Characterization Sweep")
fig.tight_layout()
plt.show()
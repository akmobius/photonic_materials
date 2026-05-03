import numpy as np
import matplotlib.pyplot as plt

# --- Configuration ---
wl = 0.785 # 785 nm
Lambda = 0.8 # Lattice constant (um)
n_sub = 1.45 # Substrate
n_inc = 1.0 # Air
n_Si = 3.71 # Core pillar material (Silicon at 785nm)

# --- 1. 2D Brillouin Zone & Diffraction Limits ---
def plot_brillouin_zone(ax):
    # Reciprocal lattice vectors (square lattice)
    Gx = 2 * np.pi / Lambda
    Gy = 2 * np.pi / Lambda

    # Brillouin Zone boundaries
    BZ_max = np.pi / Lambda
    
    # Draw BZ
    bz_square = plt.Rectangle((-BZ_max, -BZ_max), 2*BZ_max, 2*BZ_max, 
                              fill=False, edgecolor='blue', linestyle='--', lw=2, label='1st Brillouin Zone')
    ax.add_patch(bz_square)
    
    # Draw Reciprocal Lattice Points
    for m in [-1, 0, 1]:
        for n in [-1, 0, 1]:
            ax.plot(m*Gx, n*Gy, 'ko', markersize=6)
            if m==0 and n==0:
                ax.text(m*Gx+0.5, n*Gy+0.5, r'$\Gamma$', fontsize=12)
            else:
                ax.text(m*Gx+0.5, n*Gy+0.5, f'G({m},{n})', fontsize=10)

    # Draw Light Cones
    # The light cone radius is k_0 = 2*pi/wl
    k0 = 2 * np.pi / wl
    light_circle = plt.Circle((0, 0), k0, color='red', fill=False, lw=2, label='Light Circle (Air)')
    sub_circle = plt.Circle((0, 0), k0 * n_sub, color='green', fill=False, lw=2, label='Light Circle (Substrate)')
    ax.add_patch(light_circle)
    ax.add_patch(sub_circle)
    
    # Shifted light circles for diffraction orders
    # Only showing for (1,0) and (0,1) for clarity
    for m, n in [(1,0), (-1,0), (0,1), (0,-1)]:
        if m != 0 or n != 0:
            ax.add_patch(plt.Circle((m*Gx, n*Gy), k0 * n_sub, color='green', fill=False, lw=1, linestyle=':'))
        
    ax.set_aspect('equal')
    # Set limits slightly larger than BZ
    ax.set_xlim(-2*Gx, 2*Gx)
    ax.set_ylim(-2*Gy, 2*Gy)
    ax.set_title("2D Brillouin Zone & Light Cones")
    ax.set_xlabel(r"$k_x$ ($\mu m^{-1}$)")
    ax.set_ylabel(r"$k_y$ ($\mu m^{-1}$)")
    ax.legend(loc='upper right')

# --- 2. Effective Medium Theory (EMT) vs FP n_eff ---
def plot_emt_comparison(ax):
    # Pillar widths (diameters) from 0.2*Lambda to 0.8*Lambda
    diameters = np.linspace(0.2*Lambda, 0.8*Lambda, 200)
    radii = diameters / 2
    
    # Cylindrical Fill Factor
    f = (np.pi * radii**2) / (Lambda**2)
    
    # EMT Models
    # Volume Averaging bounds
    n_TE = np.sqrt(f * n_Si**2 + (1 - f) * n_inc**2) # Upper bound
    n_TM = 1 / np.sqrt(f / n_Si**2 + (1 - f) / n_inc**2) # Lower bound
    
    # Simulated FP n_eff mapping (linearly spaced 1.2 to 3.2 for the available widths)
    n_eff_FP_assumed = 1.2 + 2.0 * (diameters - 0.2*Lambda) / (0.6*Lambda)
    
    ax.plot(diameters, n_TE, 'r-', lw=2, label='EMT Volume Avg (Upper Bound)')
    ax.plot(diameters, n_TM, 'b-', lw=2, label='EMT Volume Avg (Lower Bound)')
    ax.plot(diameters, n_eff_FP_assumed, 'k--', lw=2, label='Locally Periodic $n_{eff}$ (Assumed)')
    
    # Highlight the range we used in the FP model
    ax.fill_between(diameters, 1.2, 3.21, color='gray', alpha=0.2, label='FP Valid $n_{eff}$ Range')

    ax.set_title("Effective Medium Theory vs FP $n_{eff}$")
    ax.set_xlabel(r"Pillar Diameter ($\mu$m)")
    ax.set_ylabel("Effective Refractive Index ($n_{eff}$)")
    ax.legend()
    ax.grid(True, alpha=0.3)

# --- 3. Evanescent Wave Coupling & Decay Constant ---
def plot_evanescent_coupling(ax):
    diameters = np.linspace(0.2*Lambda, 0.8*Lambda, 200)
    
    # Gap between pillars
    gap_d = Lambda - diameters
    
    # Assuming the mode effective index scales linearly for this demonstration 
    n_eff_assumed = 1.2 + 2.0 * (diameters - 0.2*Lambda) / (0.6*Lambda)
    
    # Decay constant in air: kappa = k0 * sqrt(n_eff^2 - n_air^2)
    k0 = 2 * np.pi / wl
    # Protect against sqrt(<0) if n_eff drops below n_inc (which it doesn't here)
    kappa = k0 * np.sqrt(np.maximum(n_eff_assumed**2 - n_inc**2, 0.001))
    
    # Decay length
    decay_length = 1 / kappa
    
    ax.plot(diameters, gap_d, 'k-', lw=2, label='Gap between Pillars ($d$)')
    ax.plot(diameters, decay_length, 'm--', lw=2, label=r'Evanescent Decay Length ($1/\kappa$)')
    
    # Fill region where coupling is strong (gap < decay length)
    ax.fill_between(diameters, gap_d, decay_length, where=(gap_d < decay_length), 
                    color='red', alpha=0.3, label='Strong Evanescent Coupling')
    
    ax.set_title("Evanescent Wave Coupling Analysis")
    ax.set_xlabel(r"Pillar Diameter ($\mu$m)")
    ax.set_ylabel(r"Distance ($\mu$m)")
    ax.legend()
    ax.grid(True, alpha=0.3)

# --- Main Plotting ---
if __name__ == "__main__":
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    
    plot_brillouin_zone(axs[0])
    plot_emt_comparison(axs[1])
    plot_evanescent_coupling(axs[2])
    
    plt.tight_layout()
    plt.show()

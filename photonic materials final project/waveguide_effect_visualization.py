import numpy as np
import matplotlib.pyplot as plt

# Define a simple function to approximate fundamental mode profile E(x)
def illustrative_mode_profile(x, radius, confinement_factor):
    """
    Creates an illustrative 1D electric field mode profile.
    Inside the core: Cosine-like shape.
    Outside the core: Exponential decay.
    confinement_factor (0 to 1) dictates how fast it decays outside the pillar.
    """
    E = np.zeros_like(x)
    
    # Inside the pillar (Core)
    core_mask = np.abs(x) <= radius
    # The field curves down towards the boundary based on confinement
    core_val_at_boundary = np.cos(np.pi/2 * confinement_factor)
    E[core_mask] = np.cos((x[core_mask] / radius) * (np.pi/2 * confinement_factor))
    
    # Outside the pillar (Air Cladding)
    clad_mask = np.abs(x) > radius
    # Exponential decay matching boundary value
    decay_rate = 3.0 / (1.001 - confinement_factor) # higher confinement -> much faster decay
    E[clad_mask] = core_val_at_boundary * np.exp(-decay_rate * (np.abs(x[clad_mask]) - radius))
    
    # Normalize to peak amplitude
    E = E / np.max(E)
    return E

# Create spatial x-axis
x = np.linspace(-0.8, 0.8, 1000)
n_Si = 3.71
n_Air = 1.0

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# ==========================================
# Plot 1: NARROW PILLAR (Low Confinement)
# ==========================================
r_narrow = 0.15
conf_narrow = 0.4 # Low confinement, field spreads heavily into air
E_narrow = illustrative_mode_profile(x, r_narrow, conf_narrow)
n_profile_narrow = np.piecewise(x, [np.abs(x) <= r_narrow, np.abs(x) > r_narrow], [n_Si, n_Air])

ax1.plot(x, n_profile_narrow, 'k', lw=2, label='Refractive Index $n(x)$')
ax1.fill_between(x, n_Air, n_profile_narrow, color='gray', alpha=0.1)

# Scale electric field purely for visualization purposes
E_plot_narrow = n_Air + E_narrow * (n_Si - n_Air) * 0.9 
ax1.plot(x, E_plot_narrow, 'b-', lw=3, label='Electric Field Amplitude $|E(x)|$')
ax1.fill_between(x, n_Air, E_plot_narrow, color='blue', alpha=0.1)

# Annotate the "spill"
ax1.annotate('Significant field\nleaks into "fast" air', xy=(-0.3, 1.8), xytext=(-0.7, 2.5),
             arrowprops=dict(facecolor='black', arrowstyle='->'), fontsize=10)

ax1.set_title(fr"Narrow Pillar (Radius = {r_narrow} $\mu$m)\nLow $n_{{eff}}$ (~1.5)")
ax1.set_xlabel(r"Radial Position x ($\mu$m)")
ax1.set_ylabel("Refractive Index Scale / Field Amplitude")
ax1.set_ylim(0.8, 4.0)
ax1.legend(loc='upper right')

# ==========================================
# Plot 2: WIDE PILLAR (High Confinement)
# ==========================================
r_wide = 0.35
conf_wide = 0.85 # High confinement, field mostly contained in silicon core
E_wide = illustrative_mode_profile(x, r_wide, conf_wide)
n_profile_wide = np.piecewise(x, [np.abs(x) <= r_wide, np.abs(x) > r_wide], [n_Si, n_Air])

ax2.plot(x, n_profile_wide, 'k', lw=2, label='Refractive Index $n(x)$')
ax2.fill_between(x, n_Air, n_profile_wide, color='gray', alpha=0.1)

# Scale electric field
E_plot_wide = n_Air + E_wide * (n_Si - n_Air) * 0.9 
ax2.plot(x, E_plot_wide, 'r-', lw=3, label='Electric Field Amplitude $|E(x)|$')
ax2.fill_between(x, n_Air, E_plot_wide, color='red', alpha=0.1)

# Annotate the confinement
ax2.annotate('Field tightly confined\ninside "slow" Silicon', xy=(0.3, 2.5), xytext=(0.4, 3.2),
             arrowprops=dict(facecolor='black', arrowstyle='->'), fontsize=10)

ax2.set_title(fr"Wide Pillar (Radius = {r_wide} $\mu$m)\nHigh $n_{{eff}}$ (~3.2)")
ax2.set_xlabel(r"Radial Position x ($\mu$m)")
ax2.set_ylim(0.8, 4.0)
ax2.legend(loc='upper right')

plt.suptitle("The Waveguide Effect: How Pillar Radius Dictates the Effective Index ($n_{eff}$)", fontsize=16, fontweight='bold')
plt.tight_layout()
plt.show()

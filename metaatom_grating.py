import numpy as np
import matplotlib.pyplot as plt

def get_grating_orders(wl, spacing, n_inc=1.0, n_trans=1.45, theta_inc=0):
    """
    Calculates diffraction order angles based on the grating equation.
    Ref: Lecture 12 [Source: 6309, 7901]
    """
    orders = np.arange(-3, 4)
    sin_theta_m = (n_inc * np.sin(np.deg2rad(theta_inc)) + orders * wl / spacing) / n_trans
    
    # Valid propagating orders satisfy |sin(theta)| <= 1 [Source: 6318, 7945]
    mask = np.abs(sin_theta_m) <= 1
    valid_orders = orders[mask]
    angles = np.rad2deg(np.arcsin(sin_theta_m[mask]))
    
    return valid_orders, angles

# --- Grating to Fabry-Perot Connection ---
# The individual meta-atoms act as local Fabry-Perot cavities (see metaatom_fabryperot.py).
# Each cavity imprints a specific phase shift and transmission amplitude based on its effective index.
# When these meta-atoms are arranged periodically with spacing 'pillar_spacing', they form a grating.
# The grating spacing determines the *angles* of the propagating diffraction orders, 
# while the Fabry-Perot resonance of the individual meta-atoms determines the *efficiency* 
# (how much light couples into each specific order).

# --- Configuration ---
wl = 0.785 # 785nm
pillar_spacing = 0.8 # Spacing equivalent to pillar thickness [Requested]
n_sub = 1.45

# Calculate Orders
orders, angles = get_grating_orders(wl, pillar_spacing, n_trans=n_sub)

# --- Visualization of the Array Effects & Geometry ---
fig, (ax_geom, ax_orders) = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [1, 2]})

# 1. Geometry Plot (Side View)
focal_length = 15.0 # Illustrative focal length in um
num_pillars = 30
for i in range(-num_pillars, num_pillars + 1):
    x_center = i * pillar_spacing
    
    # Calculate required phase shift for a lens profile
    phase = (2 * np.pi / wl) * (np.sqrt(x_center**2 + focal_length**2) - focal_length)
    phase_wrapped = phase % (2 * np.pi)
    
    # Map wrapped phase to pillar width (illustrative)
    pillar_width = pillar_spacing * (0.2 + 0.6 * (phase_wrapped / (2 * np.pi)))
    
    ax_geom.add_patch(plt.Rectangle((x_center - pillar_width/2, 0), pillar_width, 0.8, color='gray', alpha=0.7))

ax_geom.axhline(0, color='black', lw=2) # Substrate boundary
ax_geom.set_xlim(-(num_pillars+1)*pillar_spacing, (num_pillars+1)*pillar_spacing)
ax_geom.set_ylim(-0.2, 1.2)
ax_geom.set_aspect('auto')
ax_geom.set_title("Grating Geometry (Side View)")
ax_geom.set_xlabel("x ($\mu$m)")
ax_geom.set_ylabel("z ($\mu$m)")

# 2. Propagating Orders Plot
ax_orders.axhline(0, color='black', lw=1)
for m, ang in zip(orders, angles):
    color = 'blue' if m == 0 else 'red'
    ax_orders.arrow(0, 0, np.sin(np.deg2rad(ang)), np.cos(np.deg2rad(ang)), 
                    head_width=0.05, color=color, label=f'Order {m}' if m<=1 else "")
    ax_orders.text(np.sin(np.deg2rad(ang))*1.1, np.cos(np.deg2rad(ang))*1.1, f'm={m}', ha='center')

ax_orders.set_title(f"Propagating Grating Orders ($\Lambda$={pillar_spacing}um, $\lambda$={wl}um)")
ax_orders.set_xlim(-1.5, 1.5); ax_orders.set_ylim(-0.1, 1.5)
ax_orders.set_aspect('equal')

plt.tight_layout()
plt.show()

print(f"Propagating orders: {orders}")

# --- FP to Grating Connection Visualization ---
try:
    from metaatom_fabryperot import get_meta_atom_results
    
    # 1. Build the "Library" from the Fabry-Perot Model
    n_eff_sweep = np.linspace(1.2, 3.2, 200)
    phases = []
    for neff in n_eff_sweep:
        _, p = get_meta_atom_results(wl, 0.8, neff, n_sub, 'TE')
        phases.append(p)
    
    phases = np.unwrap(phases)
    phases = phases - np.min(phases) # Normalize to start at 0
    # Illustrative mapping from n_eff to pillar widths
    widths_library = pillar_spacing * (0.2 + 0.6 * (n_eff_sweep - 1.2) / 2.0)
    
    # 2. Target Phase for the Grating (Metalens)
    x_positions = np.linspace(-num_pillars*pillar_spacing, num_pillars*pillar_spacing, 200)
    target_phase = (2 * np.pi / wl) * (np.sqrt(x_positions**2 + focal_length**2) - focal_length)
    target_phase_wrapped = target_phase % (2 * np.pi)
    
    # 3. Create the Visualization
    fig2, axs = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot A: Fabry-Perot Phase vs Width
    axs[0].plot(widths_library, phases, 'b', lw=2)
    axs[0].set_title("1. Fabry-Perot Unit Cell Library")
    axs[0].set_xlabel("Pillar Width ($\mu$m)")
    axs[0].set_ylabel("Transmission Phase Shift (rad)")
    axs[0].grid(True, alpha=0.3)
    
    # Plot B: Target Phase Profile
    axs[1].plot(x_positions, target_phase_wrapped, 'r', lw=2)
    axs[1].set_title("2. Target Metalens Phase Profile")
    axs[1].set_xlabel("Position x ($\mu$m)")
    axs[1].set_ylabel("Required Phase (rad, wrapped)")
    axs[1].grid(True, alpha=0.3)
    
    # Plot C: Mapped Widths
    x_discrete = np.arange(-num_pillars, num_pillars + 1) * pillar_spacing
    target_phase_discrete = (2 * np.pi / wl) * (np.sqrt(x_discrete**2 + focal_length**2) - focal_length)
    target_phase_wrapped_discrete = target_phase_discrete % (2 * np.pi)
    widths_discrete = pillar_spacing * (0.2 + 0.6 * (target_phase_wrapped_discrete / (2 * np.pi)))
    
    axs[2].stem(x_discrete, widths_discrete, basefmt=" ")
    axs[2].set_title("3. Discretized Pillar Widths")
    axs[2].set_xlabel("Position x ($\mu$m)")
    axs[2].set_ylabel("Assigned Pillar Width ($\mu$m)")
    axs[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

except ImportError:
    print("Could not import get_meta_atom_results from metaatom_fabryperot. Make sure it is in the same directory.")
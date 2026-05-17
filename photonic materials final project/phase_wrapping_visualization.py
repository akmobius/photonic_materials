import numpy as np
import matplotlib.pyplot as plt

# Spatial coordinates across the lens
x = np.linspace(0, 4.5, 2000)

# 1. Target Continuous Phase (A quadratic ramp mimicking a focusing lens)
target_phase = 0.5 * x**2 * np.pi 

# 2. Perfect 2-Pi Wrapping
wrapped_2pi = target_phase % (2 * np.pi)

# 3. Imperfect 1.5-Pi Wrapping (Your pillar library maxes out at 1.5pi)
wrapped_15pi = target_phase % (1.5 * np.pi)

fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

# ==========================================
# Panel 1: The Phase Maps
# ==========================================
axs[0].plot(x, target_phase, 'k--', lw=2, label="Ideal Math Profile (Requires massive phase shift)")
axs[0].plot(x, wrapped_2pi, 'g-', lw=3, alpha=0.8, label="Perfect $2\pi$ Library (Fresnel Wrap)")
axs[0].plot(x, wrapped_15pi, 'r-', lw=2, label="Broken $1.5\pi$ Library")
axs[0].set_title("1. The Required Phase vs. Your Pillar Library", fontsize=14, fontweight='bold')
axs[0].set_ylabel("Phase Shift (radians)")

# Highlight the max thresholds
axs[0].axhline(2*np.pi, color='green', linestyle=':', alpha=0.5)
axs[0].text(0.1, 2*np.pi + 0.3, "$2\pi$ Max Threshold", color='green', fontweight='bold')
axs[0].axhline(1.5*np.pi, color='red', linestyle=':', alpha=0.5)
axs[0].text(0.1, 1.5*np.pi + 0.3, "$1.5\pi$ Max Threshold", color='red', fontweight='bold')

axs[0].legend(loc='upper left')
axs[0].set_ylim(0, 12)

# ==========================================
# Panel 2: The Wavefront for 2-Pi (Continuous)
# ==========================================
# Snapshot of the physical electric field wave: cos(phase)
wave_2pi = np.cos(wrapped_2pi)
axs[1].plot(x, wave_2pi, 'g-', lw=3)
axs[1].set_title("2. Physical Wavefront: Perfect $2\pi$ Coverage", fontsize=14, fontweight='bold', color='darkgreen')
axs[1].set_ylabel("Electric Field\nAmplitude")

# Annotate the reset points
jump_idx_2pi = np.where(np.abs(np.diff(wrapped_2pi)) > 5)[0]
for idx in jump_idx_2pi:
    axs[1].axvline(x[idx], color='green', linestyle=':', lw=2, alpha=0.4)
    
if len(jump_idx_2pi) > 0:
    axs[1].text(x[jump_idx_2pi[0]] + 0.1, 1.2, 
                "Phase resets to 0 here.\nBut $\cos(2\pi) = \cos(0)$.\nWave is perfectly continuous!", 
                color='darkgreen', fontweight='bold', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

axs[1].set_ylim(-1.5, 2.5)

# ==========================================
# Panel 3: The Wavefront for 1.5-Pi (Discontinuous)
# ==========================================
# Snapshot of the broken wave: cos(broken_phase)
wave_15pi = np.cos(wrapped_15pi)
axs[2].plot(x, wave_15pi, 'r-', lw=3)
axs[2].set_title("3. Physical Wavefront: Broken $1.5\pi$ Coverage", fontsize=14, fontweight='bold', color='darkred')
axs[2].set_ylabel("Electric Field\nAmplitude")
axs[2].set_xlabel("Position across the Metalens ($x$)", fontsize=12, fontweight='bold')

# Annotate the broken reset points
jump_idx_15pi = np.where(np.abs(np.diff(wrapped_15pi)) > 3)[0]
for idx in jump_idx_15pi:
    axs[2].axvline(x[idx], color='red', linestyle='--', lw=2, alpha=0.5)
    
if len(jump_idx_15pi) > 0:
    axs[2].text(x[jump_idx_15pi[0]] + 0.1, 1.2, 
                "Phase resets to 0 here.\nBut $\cos(1.5\pi) \\neq \cos(0)$!\nWAVEFRONT SHATTERS.", 
                color='darkred', fontweight='bold', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

axs[2].set_ylim(-1.5, 2.5)

plt.suptitle("The Physics of Phase Wrapping: Why Metalenses Require Full $2\pi$ Phase Coverage", fontsize=18, fontweight='bold', y=0.97)
plt.tight_layout()
plt.show()

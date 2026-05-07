import numpy as np
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

# ==========================================
# 1. FP Model (Locally Periodic Approximation)
# ==========================================
# Draw 3 identical pillars to represent the "locally periodic" assumption
ax1.add_patch(plt.Rectangle((1, 0), 2, 4, color='steelblue', alpha=0.8))
ax1.add_patch(plt.Rectangle((4, 0), 2, 4, color='steelblue', alpha=0.8))
ax1.add_patch(plt.Rectangle((7, 0), 2, 4, color='steelblue', alpha=0.8))
ax1.axhline(0, color='black', lw=3) # Substrate

# Draw Virtual Walls (Periodic Boundary Conditions)
ax1.axvline(3.5, color='gray', linestyle='--', lw=2)
ax1.axvline(6.5, color='gray', linestyle='--', lw=2)

# Draw Independent plane waves flowing straight down
for x in [2, 5, 8]:
    ax1.annotate('', xy=(x, 0.2), xytext=(x, 4.5),
                 arrowprops=dict(facecolor='blue', edgecolor='blue', width=3, headwidth=10))
    ax1.text(x, 2.5, '1D Plane Wave', color='white', ha='center', va='center', rotation=90, fontweight='bold')

# Annotations
ax1.text(5, 5.2, "Virtual Periodic Boundary Conditions", ha='center', color='gray', fontweight='bold')
ax1.text(5, -1.2, 
         "ASSUMPTION:\nEach pillar thinks it is surrounded\nby infinite copies of itself.\nLight 'stays in its lane'.", 
         ha='center', fontsize=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='blue', boxstyle='round,pad=0.5'))

ax1.set_xlim(0, 10)
ax1.set_ylim(-2.0, 6)
ax1.set_title("1. Fabry-Pérot / Locally Periodic (Fast Approximation)", fontsize=14, fontweight='bold')
ax1.axis('off')

# ==========================================
# 2. FDTD Model (Rigorous Maxwell Solver)
# ==========================================
# Draw 3 varying pillars (A real metalens has varying radii)
ax2.add_patch(plt.Rectangle((1, 0), 1.5, 4, color='steelblue', alpha=0.8))
ax2.add_patch(plt.Rectangle((4, 0), 2.5, 4, color='steelblue', alpha=0.8))
ax2.add_patch(plt.Rectangle((7.5, 0), 1.0, 4, color='steelblue', alpha=0.8))
ax2.axhline(0, color='black', lw=3) # Substrate

# Main waves flowing down
for x in [1.75, 5.25, 8.0]:
    ax2.annotate('', xy=(x, 0.2), xytext=(x, 4.5),
                 arrowprops=dict(facecolor='red', edgecolor='red', width=3, headwidth=10))

# --- Cross-talk and Scattering Phenomena ---
# 1. Evanescent coupling (horizontal arrows leaking between pillars)
ax2.annotate('', xy=(3.8, 2.5), xytext=(2.7, 2.5), arrowprops=dict(arrowstyle='<->', color='purple', lw=2.5))
ax2.text(3.25, 2.8, "Evanescent\nCoupling", color='purple', ha='center', fontsize=10, fontweight='bold')

ax2.annotate('', xy=(7.3, 1.5), xytext=(6.7, 1.5), arrowprops=dict(arrowstyle='<->', color='purple', lw=2.5))

# 2. Corner scattering (light diffracted off the sharp edges)
ax2.plot([4, 3.2, 2.8], [4, 4.8, 5.2], 'r--', lw=2)
ax2.plot([6.5, 7.2, 7.8], [4, 4.6, 5.0], 'r--', lw=2)
ax2.text(7.2, 5.2, "Edge Scattering", color='red', fontsize=10, fontweight='bold', ha='center')
ax2.text(3.2, 5.3, "Edge Scattering", color='red', fontsize=10, fontweight='bold', ha='center')

# Annotations
ax2.text(5, -1.2, 
         "REALITY:\nNeighboring pillars are different sizes.\nFields leak, couple, and scatter across\nthe entire 3D structure.", 
         ha='center', fontsize=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='red', boxstyle='round,pad=0.5'))

ax2.set_xlim(0, 10)
ax2.set_ylim(-2.0, 6)
ax2.set_title("2. FDTD / Maxwell Solver (Rigorous physics)", fontsize=14, fontweight='bold')
ax2.axis('off')

# Global Title
plt.suptitle("Physics Modeling: Analytical Approximations vs. Rigorous Numerical Solvers", fontsize=18, fontweight='bold', y=0.95)
plt.tight_layout()
plt.show()

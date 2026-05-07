import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12, 12))

x_center = 5

# Define vertical centers for each step
y1 = 10.0
y2 = 6.0
y3 = 2.0

# Helper function to dynamically draw perfectly fitted boxes
def draw_step(y_center, title, body, bg_color, border_color):
    # Draw Title floating just above the body box
    ax.text(x_center, y_center + 1.1, title, ha='center', va='bottom', fontsize=15, fontweight='bold', color=border_color)
    
    # Draw Body text with an auto-fitting bounding box
    ax.text(x_center, y_center, body, ha='center', va='center', fontsize=12, linespacing=2.0, 
            bbox=dict(boxstyle="round,pad=1.5", facecolor=bg_color, edgecolor=border_color, lw=2.5))

# =======================================================
# STEP 1: FORWARD DESIGN
# =======================================================
draw_step(y1, "STEP 1: FORWARD DESIGN (Analytical Approximations)", 
          "Tools: Fabry-Pérot Model / Transfer Matrix / Local Periodicity\n"
          "Action: Sweep pillar radius, build the phase library, and map the target lens profile.\n"
          "Output: Exact (x, y) coordinates and radii for all 10,000,000 pillars in the GDS file.\n"
          "Speed: Seconds/Minutes on a standard laptop.", 
          '#e6f2ff', '#0055cc')

# Down Arrow 1
ax.annotate('', xy=(x_center, y2 + 1.6), xytext=(x_center, y1 - 1.2),
            arrowprops=dict(facecolor='gray', edgecolor='gray', width=5, headwidth=15))

# =======================================================
# STEP 2: LOCAL VALIDATION
# =======================================================
draw_step(y2, "STEP 2: LOCAL VALIDATION (Rigorous but Small-Scale)", 
          "Tools: 3D FDTD (Maxwell Solver) / RCWA\n"
          "Action: Extract a tiny 'supercell' (e.g., 5 pillars representing a localized grating)\n"
          "from Step 1 and rigorously test it to see true deflection efficiency and cross-talk.\n"
          "Speed: Minutes/Hours on a workstation.", 
          '#e6ffe6', '#004d00')

# Down Arrow 2
ax.annotate('', xy=(x_center, y3 + 1.6), xytext=(x_center, y2 - 1.2),
            arrowprops=dict(facecolor='gray', edgecolor='gray', width=5, headwidth=15))

# =======================================================
# STEP 3: FINAL CHECK
# =======================================================
draw_step(y3, "STEP 3: FINAL CHECK (Rigorous but Scaled Down)", 
          "Tools: 3D FDTD (Full Device Simulation)\n"
          "Action: Simulate a drastically scaled-down 'micro-lens' version of the design\n"
          "(e.g., 20 μm wide instead of 1 mm) to prove a tight focal spot actually forms.\n"
          "Speed: Hours/Days on a supercomputing cluster. (Required before fabrication).", 
          '#ffe6e6', '#800000')

# Formatting
ax.set_xlim(0, 10)
ax.set_ylim(0, 12.5)
ax.axis('off')

plt.suptitle("The Standard Metalens Engineering Pipeline", fontsize=18, fontweight='bold', y=0.96)
plt.title("Bridging the gap between fast approximations and slow, rigorous physics", fontsize=13, style='italic', color='dimgray', pad=20)
plt.tight_layout()
plt.show()

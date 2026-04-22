import numpy as np
import matplotlib.pyplot as plt

# Physical Constants
c = 2.99792458e8  # Speed of light (m/s)
eps0 = 8.85418782e-12  # Vacuum permittivity (F/m)

# --- Problem 4.1: Optical Parametric Amplifier ---
def n_o_bbo(wl):
    """Sellmeier for BBO Ordinary index (wl in microns)"""
    return np.sqrt(2.7359 + 0.01878 / (wl**2 - 0.01822) - 0.01354 * wl**2)

def n_e_bbo(wl):
    """Sellmeier for BBO Extraordinary index (wl in microns)"""
    return np.sqrt(2.3753 + 0.01224 / (wl**2 - 0.01667) - 0.01516 * wl**2)

def n_e_theta(wl, theta_deg):
    """Index for extraordinary wave at angle theta in BBO (microns, deg)"""
    no = n_o_bbo(wl)
    ne = n_e_bbo(wl)
    theta_rad = np.radians(theta_deg)
    return 1 / np.sqrt((np.cos(theta_rad)/no)**2 + (np.sin(theta_rad)/ne)**2)

def get_group_index_o(wl):
    """Calculate ordinary group index n_g = n - lambda * dn/dlambda"""
    h = 1e-5
    n = n_o_bbo(wl)
    dn_dwl = (n_o_bbo(wl + h) - n_o_bbo(wl - h)) / (2 * h)
    return n - wl * dn_dwl

def get_group_index_e_theta(wl, theta_deg):
    """Calculate extraordinary group index for a fixed crystal angle theta"""
    h = 1e-5
    n = n_e_theta(wl, theta_deg)
    dn_dwl = (n_e_theta(wl + h, theta_deg) - n_e_theta(wl - h, theta_deg)) / (2 * h)
    return n - wl * dn_dwl

def find_theta_type1(lp_microns, ls_microns):
    """
    Finds theta such that n_e(lp, theta)/lp = n_o(ls)/ls + n_o(li)/li
    Satisfies delta_k = 0 for e -> o + o
    """
    li_microns = 1 / (1/lp_microns - 1/ls_microns)
    if li_microns > 3.0: return None
    
    target_k = n_o_bbo(ls_microns)/ls_microns + n_o_bbo(li_microns)/li_microns
    
    t_min, t_max = 0, 90
    for _ in range(30):
        t_mid = (t_min + t_max) / 2
        if n_e_theta(lp_microns, t_mid) / lp_microns > target_k:
            t_min = t_mid
        else:
            t_max = t_mid
    
    if abs(n_e_theta(lp_microns, t_mid)/lp_microns - target_k) > 0.01:
        return None
    return t_mid

def find_theta_type2(lp_microns, ls_microns):
    """
    Finds theta such that n_e(lp, theta)/lp = n_o(ls)/ls + n_e(li, theta)/li
    Satisfies delta_k = 0 for e -> o + e
    """
    li_microns = 1 / (1/lp_microns - 1/ls_microns)
    if li_microns > 3.0: return None
    
    k_s = n_o_bbo(ls_microns)/ls_microns
    
    t_min, t_max = 0, 90
    for _ in range(30):
        t_mid = (t_min + t_max) / 2
        k_p = n_e_theta(lp_microns, t_mid) / lp_microns
        k_i = n_e_theta(li_microns, t_mid) / li_microns
        if k_p - k_i > k_s:
            t_min = t_mid
        else:
            t_max = t_mid
            
    if abs(n_e_theta(lp_microns, t_mid)/lp_microns - (k_s + n_e_theta(li_microns, t_mid)/li_microns)) > 0.01:
        return None
    return t_mid

# --- Problem 4.1b: Generating 4 plots for 400 nm and 800 nm pumps ---

def plot_tuning_curves():
    lp_400 = 0.4
    lp_800 = 0.8
    
    # Ranges
    ls_400 = np.linspace(0.46, 1.5, 300)
    ls_800 = np.linspace(1.1, 2.5, 300)

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: 400 nm Pump - Type I
    thetas_400_t1 = [find_theta_type1(lp_400, s) for s in ls_400]
    axs[0, 0].plot(ls_400 * 1000, thetas_400_t1, 'b', linewidth=2)
    axs[0, 0].set_title(r'BBO 400 nm Pump: Type I ($e \rightarrow o + o$)')
    
    # Plot 2: 400 nm Pump - Type II
    thetas_400_t2 = [find_theta_type2(lp_400, s) for s in ls_400]
    axs[0, 1].plot(ls_400 * 1000, thetas_400_t2, 'g', linewidth=2)
    axs[0, 1].set_title(r'BBO 400 nm Pump: Type II ($e \rightarrow o + e$)')

    # Plot 3: 800 nm Pump - Type I
    thetas_800_t1 = [find_theta_type1(lp_800, s) for s in ls_800]
    axs[1, 0].plot(ls_800 * 1000, thetas_800_t1, 'r', linewidth=2)
    axs[1, 0].set_title(r'BBO 800 nm Pump: Type I ($e \rightarrow o + o$)')
    
    # Plot 4: 800 nm Pump - Type II
    thetas_800_t2 = [find_theta_type2(lp_800, s) for s in ls_800]
    axs[1, 1].plot(ls_800 * 1000, thetas_800_t2, 'm', linewidth=2)
    axs[1, 1].set_title(r'BBO 800 nm Pump: Type II ($e \rightarrow o + e$)')

    for ax in axs.flat:
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Phase Matching Angle (deg)')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(15, 60)

    plt.tight_layout()
    print("Four phase matching plots for Problem 4.1b generated.")

# Execute plotting
if __name__ == "__main__":
    plot_tuning_curves()
    plt.show()

    # --- Problem 4.1b: Generating 4 plots for 400 nm and 800 nm pumps ---

def plot_tuning_curves():
    lp_400 = 0.4
    lp_800 = 0.8
    
    ls_400 = np.linspace(0.46, 1.5, 300)
    ls_800 = np.linspace(1.1, 2.5, 300)

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))

    axs[0, 0].plot(ls_400 * 1000, [find_theta_type1(lp_400, s) for s in ls_400], 'b', linewidth=2)
    axs[0, 0].set_title(r'BBO 400 nm Pump: Type I ($e \rightarrow o + o$)')
    
    axs[0, 1].plot(ls_400 * 1000, [find_theta_type2(lp_400, s) for s in ls_400], 'g', linewidth=2)
    axs[0, 1].set_title(r'BBO 400 nm Pump: Type II ($e \rightarrow o + e$)')

    axs[1, 0].plot(ls_800 * 1000, [find_theta_type1(lp_800, s) for s in ls_800], 'r', linewidth=2)
    axs[1, 0].set_title(r'BBO 800 nm Pump: Type I ($e \rightarrow o + o$)')
    
    axs[1, 1].plot(ls_800 * 1000, [find_theta_type2(lp_800, s) for s in ls_800], 'm', linewidth=2)
    axs[1, 1].set_title(r'BBO 800 nm Pump: Type II ($e \rightarrow o + e$)')

    for ax in axs.flat:
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Phase Matching Angle (deg)')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(15, 60)

    plt.tight_layout()
    print("Four phase matching plots for Problem 4.1b generated.")

# --- Problem 4.1d: GVM Analysis ---

def plot_gvm_400nm():
    lp = 0.4 # microns
    ls_range = np.linspace(0.47, 1.5, 200)
    
    delta_sp = []
    delta_ip = []
    
    # Scale factor to convert (ng1 - ng2)/c to fs/mm
    # fs/mm = (1e15 / (c * 1000))
    scale = 1e12 / c # This converts to ps/mm. Then * 1000 for fs/mm.
    scale_fs_mm = (1e15 / (c * 1e3))

    for ls in ls_range:
        li = 1 / (1/lp - 1/ls)
        theta = find_theta_type1(lp, ls)
        
        if theta is not None:
            ng_p = get_group_index_e_theta(lp, theta)
            ng_s = get_group_index_o(ls)
            ng_i = get_group_index_o(li)
            
            delta_sp.append((ng_s - ng_p) * scale_fs_mm)
            delta_ip.append((ng_i - ng_p) * scale_fs_mm)
        else:
            delta_sp.append(np.nan)
            delta_ip.append(np.nan)

    plt.figure(figsize=(8, 6))
    plt.plot(ls_range * 1000, delta_sp, label=r'$\delta_{sp}$ (Signal-Pump)')
    plt.plot(ls_range * 1000, delta_ip, '--', label=r'$\delta_{ip}$ (Idler-Pump)')
    plt.axhline(0, color='black', linewidth=0.8)
    plt.xlabel('Signal Wavelength (nm)')
    plt.ylabel('GVM (fs/mm)')
    plt.title('Problem 4.1d: GVM for Type-I BBO (400 nm Pump)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Calculation for 700 nm
    ls_target = 0.7
    theta_target = find_theta_type1(lp, ls_target)
    ng_p_700 = get_group_index_e_theta(lp, theta_target)
    ng_s_700 = get_group_index_o(ls_target)
    gvm_700 = abs(ng_s_700 - ng_p_700) * scale_fs_mm
    
    pulse_duration = 100 # fs
    max_length = pulse_duration / gvm_700 # mm
    
    print(f"--- Problem 4.1d Results ---")
    print(f"Signal wavelength: 700 nm")
    print(f"Phase matching angle: {theta_target:.2f} deg")
    print(f"GVM (delta_sp): {gvm_700:.2f} fs/mm")
    print(f"Maximum crystal length for 100 fs pulse: {max_length:.3f} mm")

# Execute plotting
if __name__ == "__main__":
    plot_tuning_curves()
    plot_gvm_400nm()
    plt.show()

# --- Problem 4.2: Lithium Niobate (LiNbO3) Characterization and PPLN ---

def n_e_linbo3(wl):
    """
    Sellmeier for Congruent LiNbO3 Extraordinary index at 25 deg C.
    wl in microns. Fit from Zelmon et al. (1997).
    """
    return np.sqrt(4.5820 + 0.099169 / (wl**2 - 0.044432) - 0.02194 * wl**2)

def get_dispersion_linbo3(wl):
    """Calculates n_e, group index n_g, and GVD parameter D"""
    h = 1e-4
    n = n_e_linbo3(wl)
    dn_dwl = (n_e_linbo3(wl + h) - n_e_linbo3(wl - h)) / (2 * h)
    ng = n - wl * dn_dwl
    # GVD D = -(lambda/c)*(d2n/dl2)
    d2n_dwl2 = (n_e_linbo3(wl + h) - 2*n + n_e_linbo3(wl - h)) / (h**2)
    D = -(wl * 1e-6 / c) * (d2n_dwl2 * 1e12) # ps/(nm*km)
    return n, ng, D

def plot_4_2_a():
    """Characterize LiNbO3 at room temperature and identify zero GVD."""
    wl_range = np.linspace(0.5, 2.5, 500) # microns
    data = np.array([get_dispersion_linbo3(w) for w in wl_range])
    
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    
    # Refractive Indices
    axs[0].plot(wl_range, data[:, 0], label=r'Phase Index $n_e$', color='indigo', linewidth=2)
    axs[0].plot(wl_range, data[:, 1], '--', label=r'Group Index $n_g$', color='darkcyan', linewidth=2)
    axs[0].set_xlabel('Wavelength (microns)')
    axs[0].set_ylabel('Refractive Index')
    axs[0].set_title('Problem 4.2a: LiNbO3 Indices at 25 deg C')
    axs[0].legend()
    axs[0].grid(True, alpha=0.3)
    
    # Dispersion D
    axs[1].plot(wl_range, data[:, 2], color='purple', linewidth=2)
    axs[1].axhline(0, color='black', linewidth=1)
    idx_zero = np.argmin(np.abs(data[:, 2]))
    axs[1].annotate(f'Zero GVD ~{wl_range[idx_zero]:.3f} microns', 
                    xy=(wl_range[idx_zero], 0), xytext=(1.4, 200),
                    arrowprops=dict(facecolor='black', shrink=0.05))
    axs[1].set_xlabel('Wavelength (microns)')
    axs[1].set_ylabel(r'Dispersion D [ps/(nm$\cdot$km)]')
    axs[1].set_title('Group Velocity Dispersion in LiNbO3')
    axs[1].grid(True, alpha=0.3)
    plt.tight_layout()

def plot_4_2_b():
    """Compare Bandwidth of Uniform PPLN vs Chirped Structure (Wang et al. model)."""
    L = 10e-3 # 10 mm length
    lambda_p_center = 1.55 # pump center in microns
    
    # Phase mismatch delta_k = k(2w) - 2k(w)
    def calc_dk(wl_p):
        wl_sh = wl_p / 2.0
        return (2 * np.pi / wl_sh) * n_e_linbo3(wl_sh) - 2 * (2 * np.pi / wl_p) * n_e_linbo3(wl_p)

    wl_pump_range = np.linspace(1.54, 1.56, 1000)
    dk_vals = np.array([calc_dk(w) for w in wl_pump_range])
    
    # Uniform PPLN centered at 1.55 microns
    G_uniform = calc_dk(1.55)
    eff_uniform = np.sinc((dk_vals - G_uniform) * L / (2 * np.pi))**2
    
    # Linearly Chirped PPLN (CPPLN) following Wang et al. 2024
    # Spreads phase matching over a chirp range delta_G
    delta_G = 1500 # rad/m chirp
    # Simple model for CPPLN spectrum: integration of local phase matching
    eff_chirp = np.zeros_like(eff_uniform)
    for i, dk in enumerate(dk_vals):
        # The response is roughly a rect function over the chirp range
        if abs(dk - G_uniform) < delta_G / 2:
            eff_chirp[i] = (np.pi / delta_G) * (1/L) * 1e4 # Scaled for visibility
            
    plt.figure(figsize=(10, 6))
    plt.plot(wl_pump_range * 1000, eff_uniform, label='Uniform PPLN (Narrow Bandwidth)', color='blue')
    plt.plot(wl_pump_range * 1000, eff_chirp, label='Chirped PPLN (Wang et al. Broadband)', color='red', linewidth=3)
    plt.xlabel('Pump Wavelength (nm)')
    plt.ylabel('Normalized Conversion Efficiency')
    plt.title('Problem 4.2b: PPLN vs. Chirped Waveguide Bandwidth')
    plt.legend()
    plt.grid(True, alpha=0.3)

# --- Execute Analysis ---

if __name__ == "__main__":
    # Part A: Material Characterization
    plot_4_2_a()
    
    # Part B: Bandwidth and Wang et al. comparison
    plot_4_2_b()
    
    plt.show()
    print("Problem 4.2 characterization and comparison plots generated.")
# --- Problem 4.3: Two-Photon Microscopy Analysis ---

def problem_4_3_calculations():
    wavelength = 0.8 # microns
    delta_z_target = 5.0 # microns
    ratio = 0.8 # 80% concentration
    
    # Calculate required confocal parameter b
    # arctan(delta_z / 2b) = ratio * pi / 2
    b_req = (delta_z_target / 2.0) / np.tan(ratio * np.pi / 2.0)
    
    # Calculate w0: b = pi * w0^2 / lambda
    w0_req = np.sqrt(b_req * wavelength / np.pi)
    
    # Calculate NA approx (for air): NA = lambda / (pi * w0)
    na_req = wavelength / (np.pi * w0_req)
    
    print(f"--- Problem 4.3(c) Results ---")
    print(f"Target: {ratio*100}% signal in {delta_z_target} microns section")
    print(f"Required b: {b_req:.4f} microns")
    print(f"Max focal spot size (w0): {w0_req:.4f} microns")
    print(f"Required Numerical Aperture (NA): {na_req:.3f}")
    
    # Plot axial signal profile
    z_b = np.linspace(-5, 5, 400)
    s2p = 1 / (1 + z_b**2)
    
    plt.figure(figsize=(8, 5))
    plt.plot(z_b, s2p, color='green', linewidth=2)
    plt.fill_between(z_b, s2p, alpha=0.2, color='green', label='Fluorescence Signal')
    plt.axhline(0.5, color='red', linestyle='--', label='FWHM (z = ±b)')
    plt.xlabel('Normalized Axial distance (z/b)')
    plt.ylabel('Two-Photon Signal S2p (a.u.)')
    plt.title('Problem 4.3: Two-Photon Fluorescence Axial Profile')
    plt.legend()
    plt.grid(True, alpha=0.3)

# --- Main Plotting and Analysis ---
# --- Problem 4.1b: Generating curves for 400 nm and 800 nm pumps ---

# 4.2(a) LiNbO3 Plots
wl_range = np.linspace(0.5, 2.0, 500) # 500 nm to 2000 nm
ne_vals = n_e_linbo3(wl_range)
ng_vals = get_group_index_linbo3(wl_range)
gvd_vals = [get_gvd_linbo3(w) for w in wl_range]

plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(wl_range, ne_vals, label=r'Phase Index $n_e$')
plt.plot(wl_range, ng_vals, label=r'Group Index $n_g$', linestyle='--')
plt.xlabel('Wavelength (microns)')
plt.ylabel('Refractive Index')
plt.title('LiNbO3 Index at 25 deg C')
plt.legend()
plt.grid(True)

plt.figure()
plt.plot(wl_range, gvd_vals, color='purple')
plt.axhline(0, color='black', linewidth=1)
# Zero GVD Identification
idx_zero = np.argmin(np.abs(gvd_vals))
plt.annotate(f'Zero GVD ~{wl_range[idx_zero]:.3f} microns', 
             xy=(wl_range[idx_zero], 0), xytext=(1.2, 200),
             arrowprops=dict(facecolor='black', shrink=0.05))
plt.xlabel('Wavelength (microns)')
plt.ylabel('Dispersion D [ps/(nm*km)]')
plt.title('Group Velocity Dispersion in LiNbO3')
plt.grid(True)

# 4.2(b) PPLN SHG Efficiency
def shg_efficiency(L, delta_k):
    """Normalized SHG efficiency based on sinc squared mismatch"""
    phase = delta_k * L / 2
    # np.sinc(x) is sin(pi*x)/(pi*x)
    return np.sinc(phase / np.pi)**2

L_crystal = 10e-3 # 10 mm crystal
dk_range = np.linspace(-1000, 1000, 400)
eff = shg_efficiency(L_crystal, dk_range)

plt.figure()
plt.plot(dk_range, eff)
plt.xlabel(r'Phase Mismatch $\Delta k$ (rad/m)')
plt.ylabel('Normalized Efficiency')
plt.title('SHG Efficiency vs Phase Mismatch (L = 10 mm)')
plt.grid(True)

# Run 4.3 Analysis
problem_4_3_calculations()

plt.show()



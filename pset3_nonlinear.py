#problem 3
import numpy as np
import matplotlib.pyplot as plt

def get_ne_ppln(wl, T):
    # T in Celsius, wl in microns
    # F is the temperature factor from the provided equation
    F = (T - 24.5) * (T + 570.82)
    
    # Sellmeier coefficients from the image
    term1 = 5.35583 + 4.629e-7 * F
    term2 = (0.100473 + 3.862e-8 * F) / (wl**2 - (0.20692 - 0.89e-8 * F)**2)
    term3 = (100 + 2.657e-5 * F) / (wl**2 - 11.34927**2)
    term4 = -1.5334e-2 * wl**2
    
    return np.sqrt(term1 + term2 + term3 + term4)

# Wavelengths in microns
l1 = 1.064
l2 = 1.55
l3 = 1 / (1/l1 - 1/l2) # ~3.398 microns

temperatures = np.linspace(20, 250, 231)
periods = []

for T in temperatures:
    n1 = get_ne_ppln(l1, T)
    n2 = get_ne_ppln(l2, T)
    n3 = get_ne_ppln(l3, T)
    
    # Delta k calculation for DFG: k1 - k2 - k3
    delta_k = 2 * np.pi * (n1/l1 - n2/l2 - n3/l3)
    
    # First order QPM period: Lambda = 2*pi / delta_k
    periods.append(2 * np.pi / delta_k)

plt.figure(figsize=(10, 6))
plt.plot(temperatures, periods, color='blue', linewidth=2)
plt.title('PPLN Poling Period vs. Temperature for DFG (1064nm + 1550nm)')
plt.xlabel('Temperature (°C)')
plt.ylabel('Poling Period $\Lambda$ ($\mu$m)')
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.show()

import numpy as np
import matplotlib.pyplot as plt

def get_n_o(wl):
    # Sellmeier coefficients for KDP (wl in microns)
    return np.sqrt(1 + 2.259276*wl**2/(wl**2 - 0.01008956) + 13.00522*wl**2/(wl**2 - 400))

def get_n_e(wl):
    # Sellmeier coefficients for KDP (wl in microns)
    return np.sqrt(1 + 2.132668*wl**2/(wl**2 - 0.008637494) + 3.227381*wl**2/(wl**2 - 400))

# Wavelength range from 500 nm to 1200 nm
wavelengths = np.linspace(0.5, 1.2, 500)
theta_m_deg = []

for wl in wavelengths:
    no_w = get_n_o(wl)
    no_2w = get_n_o(wl/2)
    ne_2w = get_n_e(wl/2)
    
    # Calculate sin^2(theta_m)
    sin2_theta = (no_w**-2 - no_2w**-2) / (ne_2w**-2 - no_2w**-2)
    
    if 0 <= sin2_theta <= 1:
        theta_m_deg.append(np.degrees(np.arcsin(np.sqrt(sin2_theta))))
    else:
        theta_m_deg.append(np.nan)

plt.figure(figsize=(10, 6))
plt.plot(wavelengths * 1000, theta_m_deg, label='Type I SHG (o + o -> e)')
plt.title('KDP Phase Matching Angle vs. Fundamental Wavelength')
plt.xlabel('Fundamental Wavelength (nm)')
plt.ylabel('Phase Matching Angle (deg)')
plt.grid(True)
plt.legend()
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# KDP Sellmeier coefficients (wavelength in microns)
def n_o(wl):
    return np.sqrt(1 + 2.259276*wl**2/(wl**2 - 0.01008956) + 13.00522*wl**2/(wl**2 - 400))

def n_e(wl):
    return np.sqrt(1 + 2.132668*wl**2/(wl**2 - 0.008637494) + 3.227381*wl**2/(wl**2 - 400))

# Constants
d36 = 0.39  # pm/V for KDP
phi = np.radians(45)  # Maximizing orientation
wavelengths = np.linspace(0.5, 1.2, 500)  # 500nm to 1200nm
deff_list = []

for wl in wavelengths:
    nw = n_o(wl)
    no2w = n_o(wl/2)
    ne2w = n_e(wl/2)
    
    # Calculate phase matching angle theta_m
    arg = (nw**-2 - no2w**-2) / (ne2w**-2 - no2w**-2)
    if 0 <= arg <= 1:
        theta_m = np.arcsin(np.sqrt(arg))
        # Type I d_eff formula for KDP
        deff = -d36 * np.sin(theta_m) * np.sin(2*phi)
        deff_list.append(abs(deff))
    else:
        deff_list.append(np.nan)

plt.figure(figsize=(10, 6))
plt.plot(wavelengths * 1000, deff_list, color='crimson', lw=2)
plt.title('Effective Nonlinear Coefficient $d_{eff}$ for KDP (Type I SHG)')
plt.xlabel('Fundamental Wavelength (nm)')
plt.ylabel('$|d_{eff}|$ (pm/V)')
plt.grid(True, linestyle='--')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# KDP Sellmeier equations (wavelength in microns)
def n_o(wl):
    return np.sqrt(1 + 2.259276*wl**2/(wl**2 - 0.01008956) + 13.00522*wl**2/(wl**2 - 400))

def n_e(wl):
    return np.sqrt(1 + 2.132668*wl**2/(wl**2 - 0.008637494) + 3.227381*wl**2/(wl**2 - 400))

# Wavelength range (500nm to 1200nm)
wavelengths = np.linspace(0.5, 1.2, 500)
walkoff_deg = []

for wl in wavelengths:
    # Get indices for fundamental (w) and second harmonic (2w)
    nw_o = n_o(wl)
    n2w_o = n_o(wl/2)
    n2w_e = n_e(wl/2)
    
    # Calculate Phase Matching Angle (Type I: o + o -> e)
    # n_o(w) = n_e(2w, theta_m)
    arg = (nw_o**-2 - n2w_o**-2) / (n2w_e**-2 - n2w_o**-2)
    
    if 0 <= arg <= 1:
        theta_m = np.arcsin(np.sqrt(arg))
        
        # Calculate walk-off angle rho for the SHG wave (extraordinary at 2w)
        # tan(rho) = (n_e(2w,theta)^2 / 2) * [1/n_o(2w)^2 - 1/n_e(2w)^2] * sin(2*theta)
        # Since n_e(2w, theta_m) = n_o(w) at phase match:
        tan_rho = (nw_o**2 / 2) * ( (n2w_o**-2) - (n2w_e**-2) ) * np.sin(2 * theta_m)
        walkoff_deg.append(np.degrees(np.arctan(tan_rho)))
    else:
        walkoff_deg.append(np.nan)

plt.figure(figsize=(10, 6))
plt.plot(wavelengths * 1000, walkoff_deg, color='teal', linewidth=2)
plt.title('KDP Spatial Walk-off Angle ($\\rho$) vs. Fundamental Wavelength')
plt.xlabel('Fundamental Wavelength (nm)')
plt.ylabel('Walk-off Angle $\\rho$ (degrees)')
plt.grid(True, alpha=0.3)
plt.show()
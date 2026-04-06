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
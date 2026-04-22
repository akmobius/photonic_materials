import numpy as np
import matplotlib.pyplot as plt

# --- 1. Physics Parameters (Define these globally) ---
beta2 = 23.0    # Dispersion (fs^2/mm) - 
gamma = 5.0     # Kerr effect coefficient (W^-1 km^-1)
dt = 0.01       # Time step (ps)
dz = 0.01       # Spatial step (m)
G0 = 3.0        # Small-signal gain

# --- 2. Simulation Settings ---
num_round_trips = 100
L_fiber = 10.0  # Length (m)

# --- 3. NALM Parameters (for 5.4) ---
params = {
    'r': 0.5,           # Coupler ratio
    'L2': 0.5,          # Passive fiber length (m)
    'La': 0.6,          # Gain fiber length (m)
    'gain': G0,         # Gain
    'phi0': np.pi/2,    # Phase bias
    'gdd': 0.01         # Total dispersion compensation
}


# --- Physics ---
def propagate_nlse(u, dz, dt, beta2, gamma, g_norm=0):
    # Split-Step: 1/2 Dispersion -> Full Nonlinearity/Gain -> 1/2 Dispersion
    n = len(u)
    freq = np.fft.fftfreq(n, d=dt)
    omega = 2 * np.pi * freq
    
    # 1. Half-step Dispersion
    u_omega = np.fft.fft(u)
    u_omega *= np.exp(-1j * (beta2 / 2) * (omega**2) * (dz / 2))
    u = np.fft.ifft(u_omega)
    
    # 2. Nonlinearity and Gain
    u *= np.exp(1j * gamma * np.abs(u)**2 * dz + g_norm * dz)
    
    # 3. Half-step Dispersion
    u_omega = np.fft.fft(u)
    u_omega *= np.exp(-1j * (beta2 / 2) * (omega**2) * (dz / 2))
    u = np.fft.ifft(u_omega)
    return u
# --- Setup for 5.3a ---
dt = 0.1
tau = np.arange(-20, 20, dt)
N = 1 # Change to 2 for the second part of the check
u = N * (1.0 / np.cosh(tau))

# Propagation loop
dz = 0.01
z_max = np.pi
steps = int(z_max / dz)
u_evol = u.copy()

for _ in range(steps):
    u_evol = propagate_nlse(u_evol, dz, dt, beta2=-1, gamma=1, g_norm=0)

# --- Visualization for 5.3a ---
plt.figure(figsize=(6, 4))
plt.plot(tau, np.abs(u)**2, '--', label="Input")
plt.plot(tau, np.abs(u_evol)**2, label="Output (at z=pi)")
plt.title(f"Soliton Propagation N={N}")
plt.legend()
plt.show()

# --- Diagnostics ---
def get_fwhm(u, dt):
    intensity = np.abs(u)**2
    peak = np.max(intensity)
    half_max_indices = np.where(intensity >= peak / 2)[0]
    return (half_max_indices[-1] - half_max_indices[0]) * dt

def plot_results(u, dt, title_str):
    intensity = np.abs(u)**2
    # Shift to center zero frequency
    u_omega = np.fft.fftshift(np.fft.fft(u))
    psd = np.abs(u_omega)**2
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3))
    ax1.plot(intensity); ax1.set_title(f"{title_str} - Intensity")
    ax2.plot(psd); ax2.set_title(f"{title_str} - PSD")
    plt.show()

# --- Problem 5.3b Passive Fiber Check ---
# 1. Initialize
dt = 0.01 # ps
t = np.arange(-100, 100, dt)
# --- Initialization ---
dt = 0.01 # ps
t = np.arange(-100, 100, dt)
# Energy normalization
# Target energy E_target = 1e-9 J (1 nJ), Gaussian: u(t) = A0 * exp(-(t^2)/(2*sigma^2)),nergy of Gaussian = A0^2 * sigma * sqrt(pi)
sigma = 1.0 # Standard deviation for 1ps intensity FWHM
# FWHM of Gaussian intensity is 2 * sqrt(ln(2)) * sigma
# To get 1ps intensity FWHM, use sigma = 1 / (2 * sqrt(ln(2))) approx 0.424
sigma = 1.0 / (2 * np.sqrt(np.log(2)))
E_target = 1e-9 # 1 nJ
A0 = np.sqrt(E_target / (sigma * np.sqrt(np.pi)))
# Initialize pulse
u = A0 * np.exp(-(t**2) / (2 * sigma**2))
#u = np.exp(-(t**2) / (2 * (1.0**2))) # 1ps Gaussian
# 2. Propagate
dz = 0.01; L = 10.0 # 10 meters
for _ in range(int(L/dz)):
    u = propagate_nlse(u, dz, dt, beta2=23, gamma=5, g_norm=0)

# 3. Output
duration = get_fwhm(u, dt)
print(f"Problem 5.3b Output Duration: {duration:.4f} ps")
plot_results(u, dt, "Problem 5.3b")

# --- Problem 5.4 NALM Logic (Structural Skeleton) ---
def run_nalm_roundtrip(u, params):
    # Split
    u_cw = 1j * np.sqrt(0.5) * u
    u_ccw = -np.sqrt(0.5) * u
    
    # Propagate (CW: L2 then La; CCW: La then L2)
    u_cw = propagate_nlse(u_cw, params['L2'], dt, beta2, gamma, g_norm=0)
    u_cw = propagate_nlse(u_cw, params['La'], dt, beta2, gamma, g_norm=params['gain'])
    
    u_ccw = propagate_nlse(u_ccw, params['La'], dt, beta2, gamma, g_norm=params['gain'])
    u_ccw = propagate_nlse(u_ccw, params['L2'], dt, beta2, gamma, g_norm=0)
    
    # Recombine
    return 1j * np.sqrt(0.5) * u_cw - np.sqrt(0.5) * u_ccw


# --- Evolution Storage for 5.4 b/d ---

# 1. INITIALIZATION: Define the pulse 'u' and parameters here
dt = 0.01
t = np.arange(-100, 100, dt)
sigma = 1.0 / (2 * np.sqrt(np.log(2)))
E_target = 1e-9 
A0 = np.sqrt(E_target / (sigma * np.sqrt(np.pi)))
u = A0 * np.exp(-(t**2) / (2 * sigma**2)) # 'u' is now defined globally

# 2. PROPAGATION:
dt = 0.01
beta2 = 23.0  
gamma = 5.0   
def run_simulation(u_input):
    dz = 0.01 
    return propagate_nlse(u_input, dz, dt, beta2, gamma, g_norm=0)

u_final = run_simulation(u)

num_round_trips = 150
evolution_data = [] # To store |u|^2 for each round trip
# NALM Loop
for m in range(num_round_trips):
    u = run_nalm_roundtrip(u, params) 
    evolution_data.append(np.abs(u)**2)


# --- Heatmap for Time Evolution ---
plt.figure(figsize=(8, 6))
# Plotting the evolution of intensity over round trips
plt.imshow(np.array(evolution_data), aspect='auto', extent=[t[0], t[-1], 0, num_round_trips])
plt.colorbar(label='Normalized Power |u|^2')
plt.title("Time Evolution of Pulse over 150 Round Trips")
plt.xlabel("Time (ps)")
plt.ylabel("Round Trip Number")
plt.show()
# --- Final State Plotting ---
final_u = np.sqrt(evolution_data[-1]) # Reconstruct field from intensity
# P = |u|^2 * (scale_factor) where scale_factor comes from  normalization (let's say 1.5)
final_P = np.abs(final_u)**2 * (1.5) 

plt.figure(figsize=(8, 4))
plt.plot(t, final_P)
plt.title("Output Physical Power P(t) [W]")
plt.xlabel("Time (ps)"); plt.ylabel("Power (W)")
plt.show()

# --- Tracking Metrics ---
durations = []
peak_powers = []

for m in range(num_round_trips):
    u = run_nalm_roundtrip(u, params, dt, dz)
    durations.append(get_fwhm(u, dt))
    peak_powers.append(np.max(np.abs(u)**2))

# --- Plotting Convergence ---
plt.figure(figsize=(8, 4))
plt.plot(durations)
plt.title("Convergence of Pulse Duration")
plt.xlabel("Round Trip Number")
plt.ylabel("Pulse Duration (ps)")
plt.show()

# 1. PARAMETERS
beta2 = 23e-3  # fs^2/um (Adjusted for mm/m consistency)
gamma = 5e-3   # W^-1 m^-1
G0 = 3.0       # Gain
E_sat = 1e-9   # 1 nJ saturation energy

# 2. UPDATED NALM ROUND TRIP
def run_nalm_roundtrip(u, params, dt, dz):
    # Coupler (r=0.5)
    u_cw = 1j * np.sqrt(0.5) * u
    u_ccw = -np.sqrt(0.5) * u
    
    # Gain update
    energy = np.sum(np.abs(u)**2) * dt
    current_gain = params['gain'] / (1 + energy / E_sat)
    
    # CW: Passive L2 -> Gain La
    u_cw = propagate_nlse(u_cw, params['L2'], dt, beta2, gamma, g_norm=0)
    u_cw = propagate_nlse(u_cw, params['La'], dt, beta2, gamma, g_norm=current_gain)
    
    # CCW: Gain La -> Passive L2
    u_ccw = propagate_nlse(u_ccw, params['La'], dt, beta2, gamma, g_norm=current_gain)
    u_ccw = propagate_nlse(u_ccw, params['L2'], dt, beta2, gamma, g_norm=0)
    
    # Bias (Crucial: pi/2 enables switching)
    u_ccw *= np.exp(1j * params['phi0'])
    
    return 1j * np.sqrt(0.5) * u_cw - np.sqrt(0.5) * u_ccw

# 3. EXECUTION AND DATA COLLECTION
evolution_data = []
params = {'L2': 0.5, 'La': 0.6, 'gain': G0, 'phi0': np.pi/2}

for m in range(100):
    u = run_nalm_roundtrip(u, params, dt, dz)
    evolution_data.append(np.abs(u)**2)

# 4. PLOTTING (Part b & d)
final_u = np.sqrt(evolution_data[-1])

# Plot b: Evolution Heatmap
plt.figure(figsize=(6,4))
plt.imshow(np.array(evolution_data), aspect='auto', extent=[t[0], t[-1], 0, 100])
plt.colorbar(label='Power |u|^2')
plt.title("Evolution: Round Trips 0 to 100")
plt.show()

# Plot d: Physical Output
plt.figure(figsize=(6,4))
plt.plot(t, np.abs(final_u)**2) # Ensure this is not 0!
plt.title("Final Output Pulse P(t)")
plt.xlabel("Time (ps)"); plt.ylabel("Power (W)")
plt.show()

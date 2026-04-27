import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def ricker_wavelet(t, f0=2.0, t0=0.5):
    """
    Ricker wavelet source time function (Force magnitude).
    f0: Dominant frequency
    t0: Time shift (to center the peak)
    """
    arg = (np.pi * f0 * (t - t0))**2
    return (1.0 - 2.0 * arg) * np.exp(-arg)

def calculate_displacement(x, y, t, rho, cp, cs, f0, t0, force_dir=(1, 0)):
    """
    Calculates the displacement field u(x,t) based on the Stokes solution.
    
    Parameters:
    x, y : 2D meshgrid arrays of coordinates
    t    : Current simulation time
    rho  : Density
    cp   : P-wave velocity
    cs   : S-wave velocity
    force_dir : Tuple (fx, fy) direction of the point force
    """
    
    # 1. Geometry and Distance
    r = np.sqrt(x**2 + y**2)
    # Avoid division by zero at the source
    epsilon = 1e-6
    r_safe = np.maximum(r, epsilon)
    
    # Direction cosines (gamma)
    gamma_x = x / r_safe
    gamma_y = y / r_safe
    
    # Force direction components (beta index in the formula)
    # We assume the force is in the x-y plane.
    # If force is pure x-direction: fx=1, fy=0.
    fx, fy = force_dir
    
    # 2. Retarded Times
    t_p = t - r_safe / cp
    t_s = t - r_safe / cs
    
    # 3. Source Function Evaluations for Far-Field
    # F(t - r/c)
    F_p = ricker_wavelet(t_p, f0, t0)
    F_s = ricker_wavelet(t_s, f0, t0)
    
    # 4. Near-Field Integral Term
    # Integral_{r/cp}^{r/cs} tau * F(t - tau) dtau
    # We compute this numerically using the trapezoidal rule.
    # Create a 3rd dimension for integration steps
    num_steps = 30
    # Shape: (num_steps, height, width)
    tau_start = r_safe / cp
    tau_end = r_safe / cs
    
    # Broadcasting tau over the grid
    # We create a normalized integration variable s in [0, 1]
    s = np.linspace(0, 1, num_steps)[:, None, None] 
    
    # Map s to tau: tau = start + s * (end - start)
    tau_grid = tau_start + s * (tau_end - tau_start)
    dt_grid = (tau_end - tau_start) / (num_steps - 1)
    
    # Evaluate integrand: tau * F(t - tau)
    F_tau = ricker_wavelet(t - tau_grid, f0, t0)
    integrand = tau_grid * F_tau
    
    # Integrate over the 0-th axis (num_steps)
    integral_term = np.trapz(integrand, axis=0, dx=1.0) * dt_grid[0] 
    # Note: dt_grid[0] is the dt for the first pixel, but dt varies by pixel.
    # Correct way for varying integration limits in numpy:
    # Since interval length (tau_end - tau_start) varies by pixel, 
    # the 'dx' in trapz should actually be the interval length / steps.
    # A simpler vectorization:
    integral_val = np.sum(integrand, axis=0) * (tau_end - tau_start) / num_steps

    # 5. Assemble the Terms
    
    # Initialize displacement components
    ux = np.zeros_like(x)
    uy = np.zeros_like(y)
    
    # Constants
    coef_p = 1.0 / (4 * np.pi * rho * cp**2)
    coef_s = 1.0 / (4 * np.pi * rho * cs**2)
    coef_nf = 1.0 / (4 * np.pi * rho)
    
    # We iterate over the force direction components (j or beta)
    # Since we are in 2D, j can be x or y.
    
    # --- Contribution from Force in X (beta=x) ---
    if fx != 0:
        # P-wave term (Far-field)
        # Term: gamma_i * gamma_x * F_p / r
        # i=x: gamma_x * gamma_x
        # i=y: gamma_y * gamma_x
        term_p_x = (gamma_x * gamma_x) * F_p / r_safe
        term_p_y = (gamma_y * gamma_x) * F_p / r_safe
        
        # S-wave term (Far-field)
        # Term: (delta_ix - gamma_i * gamma_x) * F_s / r
        # i=x: (1 - gamma_x^2)
        # i=y: (0 - gamma_y*gamma_x)
        term_s_x = (1.0 - gamma_x * gamma_x) * F_s / r_safe
        term_s_y = (0.0 - gamma_y * gamma_x) * F_s / r_safe
        
        # Near-field term
        # Term: (3 * gamma_i * gamma_x - delta_ix) * integral / r^3
        term_nf_x = (3 * gamma_x * gamma_x - 1.0) * integral_val / (r_safe**3)
        term_nf_y = (3 * gamma_y * gamma_x - 0.0) * integral_val / (r_safe**3)
        
        ux += fx * (coef_p * term_p_x + coef_s * term_s_x + coef_nf * term_nf_x)
        uy += fx * (coef_p * term_p_y + coef_s * term_s_y + coef_nf * term_nf_y)

    # --- Contribution from Force in Y (beta=y) ---
    if fy != 0:
        # P-wave
        term_p_x = (gamma_x * gamma_y) * F_p / r_safe
        term_p_y = (gamma_y * gamma_y) * F_p / r_safe
        
        # S-wave
        term_s_x = (0.0 - gamma_x * gamma_y) * F_s / r_safe
        term_s_y = (1.0 - gamma_y * gamma_y) * F_s / r_safe
        
        # Near-field
        term_nf_x = (3 * gamma_x * gamma_y - 0.0) * integral_val / (r_safe**3)
        term_nf_y = (3 * gamma_y * gamma_y - 1.0) * integral_val / (r_safe**3)
        
        ux += fy * (coef_p * term_p_x + coef_s * term_s_x + coef_nf * term_nf_x)
        uy += fy * (coef_p * term_p_y + coef_s * term_s_y + coef_nf * term_nf_y)

    return ux, uy

# --- Simulation Parameters ---
# Material properties (similar to the provided PDF context)
rho = 2670.0  # kg/m^3
cp = 6000.0   # m/s
cs = 3464.0   # m/s

# Source
f0 = 15.0     # Frequency (Hz)
t0 = 0.15      # Source delay (s)
force_vector = (1.0, 0.0) # Point force in x-direction

# Grid
L = 4000.0    # Domain size (meters)
N = 300       # Resolution
x = np.linspace(-L, L, N)
y = np.linspace(-L, L, N)
X, Y = np.meshgrid(x, y)

# Time snapshot
current_time = 0.6  # Seconds

# --- Compute ---
Ux, Uy = calculate_displacement(X, Y, current_time, rho, cp, cs, f0, t0, force_dir=force_vector)
Magnitude = np.sqrt(Ux**2 + Uy**2)

# --- Plot ---
plt.figure(figsize=(10, 8))
plt.imshow(Magnitude, extent=[-L, L, -L, L], origin='lower', cmap='inferno', norm=Normalize(vmin=0, vmax=np.percentile(Magnitude, 99)))
plt.colorbar(label='Displacement Magnitude |u|')
plt.title(f'Analytical Stokes Solution (t={current_time}s)\nPoint Force in X-direction')
plt.xlabel('x (m)')
plt.ylabel('y (m)')

# Annotations
plt.text(0, L*0.85, 'P-Wave', color='cyan', ha='center', fontsize=12, fontweight='bold')
plt.text(0, L*0.45, 'S-Wave', color='white', ha='center', fontsize=12, fontweight='bold')
plt.grid(True, alpha=0.3)

plt.show()
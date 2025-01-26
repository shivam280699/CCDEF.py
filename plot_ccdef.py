import matplotlib.pyplot as plt
from ccdef import CCDEF

def run_calculation():
    # Initialize CCDEF and parameters
    cc = CCDEF()
    init_params = {
        'wma': 19.0,    # Projectile mass (negative if deformed)
        'wza': 9.0,     # Projectile charge
        'wmb': 150.0,   # Target mass (negative if deformed)
        'wzb': 60.0,     # Target charge
        'dv': 20.0,      # Barrier adjustment
        'fcc': 1,        # Coupling form factor type
        'be2a': 0.0,     # Projectile quadrupole deformation
        'be4a': 0.0,     # Projectile hexadecapole deformation
        'be2b': 0.0,    # Target quadrupole deformation
        'be4b': 0.0     # Target hexadecapole deformation
    }
    cc.initialize_parameters(**init_params)

    # Define channels
    surface_channels = [
        {'beta': 0.13, 'flam': 2, 'q': -0.285},  # Projectile mode
        {'beta': 0.381, 'flam': 4, 'q': -0.05}     # Target mode
    ]
    additional_channels = []

    # Calculate cross-sections
    energies, coupled, uncoupled, _, _ = cc.calculate_cross_sections(
        emin=54.0, emax=85.0, de=2.0, ns=2, na=0,
        surface_channels=surface_channels, additional_channels=additional_channels
    )
    return energies, coupled, uncoupled

# --- Generate data and plot ---
energies, coupled, uncoupled = run_calculation()

plt.figure(figsize=(8, 6))
plt.semilogy(energies, coupled, 'b--', linewidth=2, label='Coupled Channels')
plt.semilogy(energies, uncoupled, 'r-', linewidth=2, label='Uncoupled Channels')

plt.xlabel('Energy (MeV)', fontsize=12)
plt.ylabel('Cross Section (mb)', fontsize=12)
plt.title('Fusion Cross Sections', fontsize=14)
plt.legend()
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.tight_layout()

# Save and show plot
plt.savefig('cross_sections_log.png', dpi=300)
plt.show()
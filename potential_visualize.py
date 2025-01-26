import numpy as np
import matplotlib.pyplot as plt

class CCDEF:
    def __init__(self):
        self.PI = np.pi
        self.AU = 931.5016        # Atomic mass unit (MeV/c²)
        self.HC = 197.3286        # ℏc (MeV·fm)
        self.FPI = 3.544908       # Scaling factor from Fortran
        self.idefa = 0            # Deformation flags
        self.idefb = 0

    def initialize_parameters(self, wma, wza, wmb, wzb, dv, fcc, be2a, be4a, be2b, be4b):
        self.idefa = 1 if wma < 0 else 0
        self.idefb = 1 if wmb < 0 else 0
        self.wma = abs(wma)
        self.wmb = abs(wmb)
        self.wza = wza
        self.wzb = wzb
        self.dv = dv
        self.fcc = fcc
        self.be2a = be2a if self.idefa else 0.0
        self.be4a = be4a if self.idefa else 0.0
        self.be2b = be2b if self.idefb else 0.0
        self.be4b = be4b if self.idefb else 0.0

        # Nuclear geometry calculations
        cbrt_wma = self.wma ** (1/3)
        cbrt_wmb = self.wmb ** (1/3)
        self.ra = 1.233 * cbrt_wma - 0.978 / cbrt_wma
        self.rb = 1.233 * cbrt_wmb - 0.978 / cbrt_wmb
        self.rred = (self.ra * self.rb) / (self.ra + self.rb)
        self.rab = self.ra + self.rb + 0.29
        self.redm = (self.wma * self.wmb) / (self.wma + self.wmb)
        self.h2m = (self.HC ** 2) / self.AU
        self.etak = 1.43997 * self.wza * self.wzb  # Coulomb constant
        self.a0r = 0.63

        # Nuclear potential depth
        v0rx = 30.08 * (1 - 1.8*(1 - 2*self.wza/self.wma)*(1 - 2*self.wzb/self.wmb))
        self.v0r = v0rx * self.rred + self.dv - 20.0
        self.update_deformation_params(0.0, 0.0)

    def update_deformation_params(self, theda, thedb):
        if self.idefa:
            xsqa = np.cos(theda)**2
            self.b2y2a = self.be2a * np.sqrt(5/(4*self.PI)) * 0.5*(3*xsqa - 1)
            self.b4y4a = self.be4a * np.sqrt(9/(4*self.PI)) * (35*xsqa**2 - 30*xsqa + 3)/8
        else:
            self.b2y2a = 0.0
            self.b4y4a = 0.0

        if self.idefb:
            xsqb = np.cos(thedb)**2
            self.b2y2b = self.be2b * np.sqrt(5/(4*self.PI)) * 0.5*(3*xsqb - 1)
            self.b4y4b = self.be4b * np.sqrt(9/(4*self.PI)) * (35*xsqb**2 - 30*xsqb + 3)/8
        else:
            self.b2y2b = 0.0
            self.b4y4b = 0.0

    def pot(self, rr):
        ss = rr - self.rab - self.ra*(self.b2y2a + self.b4y4a) - self.rb*(self.b2y2b + self.b4y4b)
        arg = np.exp(-ss / self.a0r)
        arg1 = 1.0 + arg
        ur = -self.v0r * arg / arg1
        dur = -ur / (self.a0r * arg1)
        ddur = dur * (1 - 2/arg1) / self.a0r
        return ur, dur, ddur

    def potent(self, rx):
        x = rx
        x2, x3, x4 = x*x, x**3, x**4
        r2a, r4a = self.ra**2, self.ra**4
        r2b, r4b = self.rb**2, self.rb**4

        v0_coul = self.etak/x * (1 + 0.6*self.b2y2a*r2a/x2 + self.b4y4a*r4a/(3*x4) +
                                  0.6*self.b2y2b*r2b/x2 + self.b4y4b*r4b/(3*x4))
        ur, dur, ddur = self.pot(x)
        v1_coul = -self.etak/x2*(1 + 1.8*self.b2y2a*r2a/x2 + 5*self.b4y4a*r4a/(3*x4) +
                                  1.8*self.b2y2b*r2b/x2 + 5*self.b4y4b*r4b/(3*x4)) + dur
        v2_coul = 2*self.etak/x3*(1 + 3.6*self.b2y2a*r2a/x2 + 15*self.b4y4a*r4a/(3*x4) +
                                   3.6*self.b2y2b*r2b/x2 + 15*self.b4y4b*r4b/(3*x4)) + ddur
        return v0_coul + ur, v1_coul, v2_coul

    def find_barrier(self):
        rmax, dr, s = 20.0, 0.5, -1.0
        rbz, y1 = rmax, -1.0
        while True:
            _, v1, v2 = self.potent(rbz)
            if v1 * y1 < 0:
                if dr < 0.01:
                    break
                dr *= 0.5
                s = -s
            rbz += s * dr
            # print(rbz)
            y1 = v1
            if rbz < 0.8 or rbz > rmax:
                return None, None, -1.0  # Invalid barrier
        v0, _, v2 = self.potent(rbz)
        homega = self.HC * np.sqrt(-v2 / (self.redm * self.AU))
        return rbz, v0, homega

    def compute_potential_profile(self, r_min=1.0, r_max=20.0, num_points=200):
        """Compute all potential components at different radii"""
        radii = np.linspace(r_min, r_max, num_points)
        v_total = []
        v_coulomb = []
        v_nuclear = []
        
        self.update_deformation_params(0.0, 0.0)  # Spherical case
        
        for r in radii:
            # Get individual components
            v0_tot, _, _ = self.potent(r)
            ur, _, _ = self.pot(r)
            v_total.append(v0_tot)
            v_coulomb.append(self.etak/r)  # Basic Coulomb term
            v_nuclear.append(ur)
            
        return radii, np.array(v_coulomb), np.array(v_nuclear), np.array(v_total) 

    def visualize_potential(self, r_min=1.0, r_max=20.0, num_points=200, save_path=None):
        """Plot all potential components with barrier parameters"""
        # Compute potentials
        radii, v_coul, v_nuc, v_tot = self.compute_potential_profile(r_min, r_max, num_points)
        
        # Find barrier parameters
        rbar, vb, homega = self.find_barrier()
        
        # Create plot
        plt.figure(figsize=(10, 7))
        
        # Plot components
        plt.plot(radii, v_coul, 'b--', lw=1.5, label='Coulomb Potential')
        plt.plot(radii, v_nuc, 'g--', lw=1.5, label='Nuclear Potential') 
        plt.plot(radii, v_tot, 'k-', lw=2, label='Total Potential')
        
        # Add barrier markers if found
        if homega > 0:
            plt.axvline(rbar, color='r', ls=':', label=f'Barrier Radius ({rbar:.2f} fm)')
            plt.axhline(vb, color='m', ls=':', label=f'Barrier Height ({vb:.2f} MeV)')
            
            # Annotate barrier point
            plt.plot(rbar, vb, 'ro', markersize=8, 
                    markeredgecolor='k', markerfacecolor='gold')
        
        # Formatting
        plt.xlabel('Radius (fm)', fontsize=12)
        plt.ylabel('Potential Energy (MeV)', fontsize=12)
        plt.title('Potential Energy Components', fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.legend(loc='upper right')
        plt.ylim(min(v_nuc)-5, max(v_coul)+5)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()

if __name__ == "__main__":
    cc = CCDEF()

    # Example input (modify as needed)
    init_params = {
        'wma': 19.0,    # Projectile mass (negative if deformed)
        'wza': 9.0,    # Projectile charge
        'wmb': 150.0,  # Target mass (negative if deformed)
        'wzb': 60.0,    # Target charge
        'dv': 20.0,     # Barrier adjustment
        'fcc': 1,       # Coupling form factor type
        'be2a': 0.0,    # Projectile quadrupole deformation
        'be4a': 0.0,    # Projectile hexadecapole deformation
        'be2b': 0.0,   # Target quadrupole deformation
        'be4b': 0.0    # Target hexadecapole deformation
    }
    cc.initialize_parameters(**init_params)

    cc.visualize_potential(r_min=8.0, r_max=20.0, save_path='potential_profile.png')
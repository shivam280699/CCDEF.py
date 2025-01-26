import numpy as np

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
            print(rbz)
            y1 = v1
            if rbz < 0.8 or rbz > rmax:
                return None, None, -1.0  # Invalid barrier
        v0, _, v2 = self.potent(rbz)
        homega = self.HC * np.sqrt(-v2 / (self.redm * self.AU))
        return rbz, v0, homega

    def calculate_cross_sections(self, emin, emax, de, ns, na, surface_channels, additional_channels):
        nmax = ns + na
        ne = int((emax - emin)/de + 1.001)
        if ne > 51:
            raise ValueError("Maximum 51 energy points allowed.")

        # Initialize result arrays
        sum_coupled = np.zeros(ne, dtype=np.float64)
        sum_uncoupled = np.zeros(ne, dtype=np.float64)
        sigl = np.zeros(100, dtype=np.float64)   # Partial wave cross-sections (coupled)
        sigl0 = np.zeros(100, dtype=np.float64)  # Partial wave cross-sections (uncoupled)

        # Angular integration parameters
        ithdma = 90 if self.idefa else 1
        ithdmb = 90 if self.idefb else 1
        pdws = 0.0

        # Angular integration loops
        for ithda in range(1, ithdma + 1):
            theda = self.PI * ithda / 180.0
            pdwa = np.sin(theda) * self.PI/180.0 if ithdma > 1 else 1.0
            
            for ithdb in range(1, ithdmb + 1):
                thedb = self.PI * ithdb / 180.0
                pdwb = np.sin(thedb) * self.PI/180.0 if ithdmb > 1 else 1.0
                pdw = pdwa * pdwb
                pdws += pdw

                # Update deformation parameters for current angles
                self.update_deformation_params(theda, thedb)
                rbar, vb, homega = self.find_barrier()
                if ns == 0 or homega < 0:
                    continue  # Skip invalid barriers

                # Process surface channels (dynamic calculation at rbar)
                channels = []
                for ch in surface_channels:
                    beta = ch['beta']
                    flam = abs(ch['flam'])
                    q = ch['q']
                    rrt = self.rb if ch['flam'] > 0 else self.ra
                    ur, dur, ddur = self.pot(rbar)
                    
                    # Coupling strength calculation
                    fs = beta * rrt * (-dur + 3*self.etak*(rrt/rbar)**(flam-1) / 
                                      (2*flam + 1) / rbar**2) / self.FPI
                    sq = np.sqrt(q**2 + 4*fs**2)
                    fla = [0.5*(-q - sq), 0.5*(-q + sq)]
                    pa = [fs**2/(fs**2 + fla[0]**2), fs**2/(fs**2 + fla[1]**2)]
                    fkap = -ddur/dur if dur != 0 else 0.0
                    dfla = [2*fkap*fs**2/sq, -2*fkap*fs**2/sq]
                    dfla2 = [
                        -4*(fkap*fs)**2/sq + 8*(fkap*fs**2)**2/sq**3,
                        4*(fkap*fs)**2/sq - 8*(fkap*fs**2)**2/sq**3
                    ]
                    channels.append({'fla': fla, 'pa': pa, 'dfla': dfla, 'dfla2': dfla2})

                # Process additional channels (if any)
                for ch in additional_channels:
                    f = ch['f']
                    q = ch['q']
                    fkap = 0.71  # Fixed value from Fortran
                    sq = np.sqrt(q**2 + 4*f**2)
                    fla = [0.5*(-q - sq), 0.5*(-q + sq)]
                    pa = [f**2/(f**2 + fla[0]**2), f**2/(f**2 + fla[1]**2)]
                    dfla = [2*fkap*f**2/sq, -2*fkap*f**2/sq]
                    dfla2 = [
                        -4*(fkap*f)**2/sq + 8*(fkap*f**2)**2/sq**3,
                        4*(fkap*f)**2/sq - 8*(fkap*f**2)**2/sq**3
                    ]
                    channels.append({'fla': fla, 'pa': pa, 'dfla': dfla, 'dfla2': dfla2})

                # Generate all channel combinations
                nd = 2 ** nmax
                npows = [2 ** (nmax - 1 - n) for n in range(nmax)]
                
                for i1 in range(nd):
                    p = 1.0
                    fl = 0.0
                    dfl = 0.0
                    dfl2 = 0.0
                    ick = 0
                    
                    for n in range(nmax):
                        n_bit = (i1 // npows[n]) % 2
                        ick += n_bit
                        if ick > nmax:
                            break
                        ch = channels[n]
                        p *= ch['pa'][n_bit]
                        fl += ch['fla'][n_bit]
                        dfl += ch['dfla'][n_bit]
                        dfl2 += ch['dfla2'][n_bit]
                    
                    if ick > nmax:
                        continue

                    # Barrier adjustment
                    delta = self.fcc * dfl / (self.redm * (homega**2)/self.h2m - dfl2)
                    if abs(delta) > 0.99:
                        delta = -0.99 * np.sign(fl) if fl != 0 else 0.0
                    vbw = vb - 0.5*self.redm*(homega*delta)**2/self.h2m + fl + dfl*delta + 0.5*dfl2*delta**2
                    rcald = rbar + delta
                    eps = homega / (2 * self.PI)
                    facwd = 31.416 * rcald**2 * eps
                    facw = 31.416 * rbar**2 * eps

                    # Energy loop
                    for ie in range(ne):
                        e = emin + de * ie
                        aqa = (e - vbw)/eps
                        if aqa < 30.0:
                            aqa = np.log1p(np.exp(aqa))
                        sum_coupled[ie] += pdw * facwd * p * aqa / e

                        # Uncoupled calculation (i1=0)
                        if i1 == 0:
                            aqa0 = (e - vb)/eps
                            if aqa0 < 30.0:
                                aqa0 = np.log1p(np.exp(aqa0))
                            sum_uncoupled[ie] += pdw * facw * aqa0 / e

                        # Partial wave calculation (only for single energy)
                        if ne == 1:
                            for il in range(100):
                                gl = il  # l = 0,1,...,99
                                vbwl = vbw + 0.5*self.HC**2*gl*(gl+1)/(self.redm*self.AU*rcald**2)
                                vbl = vb + 0.5*self.HC**2*gl*(gl+1)/(self.redm*self.AU*rbar**2)
                                factor = 31.41592*(2*gl+1)*self.HC**2/(2*self.redm*self.AU*e)
                                
                                # Transmission coefficients
                                aux = np.exp((e - vbwl)/eps)
                                tl = np.exp(aux)/(1 + np.exp(aux)) if aux < 30.0 else 1.0
                                sigl[il] += pdw * factor * p * tl
                                
                                if i1 == 0:
                                    aux0 = np.exp((e - vbl)/eps)
                                    tl0 = np.exp(aux0)/(1 + np.exp(aux0)) if aux0 < 30.0 else 1.0
                                    sigl0[il] += pdw * factor * tl0

        # Normalize results
        sum_coupled /= pdws
        sum_uncoupled /= pdws
        sigl /= pdws
        sigl0 /= pdws
        
        energies = np.linspace(emin, emax, ne)
        return energies, sum_coupled, sum_uncoupled, sigl, sigl0

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

    # Define channels (modify as needed)
    surface_channels = [
        {'beta': 0.13, 'flam': 2, 'q': -0.285},
        {'beta': 0.381, 'flam': 4, 'q': -0.05}
    ]
    additional_channels = []  # No additional channels

    # Run calculation
    energies, coupled, uncoupled, sigl, sigl0 = cc.calculate_cross_sections(
        emin=54.0, emax=85.0, de=2.0, ns=2, na=0,
        surface_channels=surface_channels, additional_channels=additional_channels
    )

    # Print cross-sections
    if len(energies) > 1:
        print("\nFusion cross sections:")
        print("Energy (MeV) | Coupled (mb) | Uncoupled (mb)")
        print("-" * 45)
        for e, c, uc in zip(energies, coupled, uncoupled):
            print(f"{e:11.1f} | {c:12.2e} | {uc:12.2e}")
    else:
        print(f"\nCross section at E = {energies[0]:.1f} MeV:")
        print(f"Coupled:   {coupled[0]:.2e} mb")
        print(f"Uncoupled: {uncoupled[0]:.2e} mb")

    # Print partial waves and moments for single energy
    if len(energies) == 1:
        print("\nPartial wave contributions (mb/ħ):")
        print(" l | Coupled        | Uncoupled")
        print("-" * 35)
        for il in range(100):
            l = il
            if sigl[il] < 1.2e-6:
                continue  # Skip small contributions
            print(f"{l:2} | {sigl[il]:14.4e} | {sigl0[il]:14.4e}")

        # Calculate moments
        s0 = np.sum(sigl)
        s1 = np.sum([l * sigl[l] for l in range(100)])
        s2 = np.sum([l**2 * sigl[l] for l in range(100)])
        su0 = np.sum(sigl0)
        su1 = np.sum([l * sigl0[l] for l in range(100)])
        su2 = np.sum([l**2 * sigl0[l] for l in range(100)])

        avl = s1 / s0 if s0 != 0 else 0.0
        sd = np.sqrt(s2/s0 - avl**2) if s0 != 0 else 0.0
        avl0 = su1 / su0 if su0 != 0 else 0.0
        sd0 = np.sqrt(su2/su0 - avl0**2) if su0 != 0 else 0.0

        print("\nMoments of the angular momentum distribution:")
        print(f"Coupled:   <l> = {avl:.2e}, σ = {sd:.2e}")
        print(f"Uncoupled: <l> = {avl0:.2e}, σ = {sd0:.2e}")

    # Write results to file (optional)
    with open("cross_sections.dat", "w") as f:
        f.write("Energy (MeV) | Coupled (mb) | Uncoupled (mb)\n")
        for e, c, uc in zip(energies, coupled, uncoupled):
            f.write(f"{e:.1f} \t\t {c:.4e} \t\t {uc:.4e}\n")
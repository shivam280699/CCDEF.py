 
## <u>Input files</u>

`wma, wmb`: mass of projectile and target respectively
`wza, wzb`: atomic number of projectile and target respectively
`dv`: parameter to adjust barrier parameters (typically -10 to +10)
`fcc`: coupling form factor type (0 or 1) 
`be2a`: static quadruple deformation for projectile
`be4a`: static hexadecapole deformation for projectile
`be2b`: static quadruple deformation for target
`be4b`: static hexadecapole deformation for target

`emin`: minimum energy for cross section
`emax`: maximum energy for cross section
`de`: energy interval

`ns`: number of surface inelastic channels
`na`: number of additional channels

## <u>Code</u>

`sigl0`: Uncoupled-channel partial wave cross-section for momentum l = il - 1
`sigl`: Coupled (bare nucleus) partial wave cross-section
`sum0`: Total Uncoupled-channel cross-section at energy e = emin + de*(ie-1)
`sum`: Total coupled cross-section

`idefa`: Flag indicating projectile deformation (1 if projectile mass is negative)
`idefb`: Flag indicating target deformation (1 if target mass is negative)

`nmax`: Maximum number of coupling channels
`ne`: Number of energy points to calculate

`ra`: Nuclear radius of the projectile
`rb`: Nuclear radius of the target

`rred`: Reduced radius (geometric mean of projectile and target radii)
`rab`: Total interaction distance
`redm`: Reduced mass of the nuclear system

`au`: atomic mass unit in MeV/$c^{2}$
`hc`: product of plank's constant and speed of light in MeV-fm
`h2m`: Reduced Planck's constant squared divided by nuclear mass (energy-mass conversion factor)
`fpi`: 
`a0`: diffuseness parameter in fm


`etak`: coulomb parameter $k_e$. 
This term calculates the Coulomb potential energy between two nuclei at a separation distance 
$V_{\text{Coulomb}} = \frac{k_e Z_1 Z_2}{r}$
($Z_1$, $Z_2$): Charges of the projectile and target.  
$k_e$ = $\frac{e^2}{4\pi\epsilon_0}$ $\approx$ 1.43997 : Fundamental constant for electrostatic interactions. 


`pdws`: probability wrights 

`ithdma`: iterator variable used for angular integration
`theda`: Angular variable for projectile deformation angle (in radians)
`pdwa`: Probability density weight for projectile angular distribution

`xsqa`: Cosine squared of the projectile deformation angle
`b2y2a`: Modified quadrupole deformation parameter for projectile
`b4y4a`: Modified hexadecapole deformation parameter for projectile

`ithdmb`: iterator variable used for target angular integration
`thedb`: Angular variable for target deformation angle (in radians)
`pdwb`: Probability density weight for target angular distribution 

`xsqb`: Cosine squared of the target deformation angle
`b2y2b`: Modified quadrupole deformation parameter for target
`b4y4b`: Modified hexadecapole deformation parameter for target

`pdw`: Combined probability density weight (projectile * target)

`v0rx`: Initial nuclear potential depth
$v_{0rx} = 30.08 \cdot \left(1 - 1.8 \cdot (1 - \frac{2Z_p}{M_p}) \cdot (1 - \frac{2Z_t}{M_t})\right)$
`v0r`: Adjusted nuclear potential depth
$v_{0r} = v_{0rx} \cdot r_{red} + dv - 20$

`eps`: Energy scale parameter
`rcal/rbar`: Characteristic nuclear interaction radius 

`rrt`: Radius used for specific nuclear channel calculations
`flam/fslam`: Multipolarity of nuclear deformation
`fs`: Coupling strength for nuclear channels
`fkap`: Potential gradient parameter

`sq`: Square root of coupling-related calculation
`fla`: Channel energy levels
`dfla`: First derivative of channel energy levels
`dfla2`: Second derivative of channel energy levels
`pa`: Probability of being in a specific channel 

`delta`: Coupling-induced energy shift
`vbw`: Modified barrier height
`rcald`: Modified interaction radius
`facw/facwd`: Scaling factors for cross-section calculations 

`aqa`: Exponential energy scaling factor

`gl`: Angular momentum quantum number
`vbwl/vbl`: Channel-specific barrier heights
`aux`: Auxiliary calculation variable
`sigl/sigl0`: Partial wave cross-sections (coupled and uncoupled) 

`pdws`: Total probability density weight

## <u>Function pot</u>

`rr`: radial distance
`ur`: nuclear potential
`dur`: first derivative of the potential
`ddur`: second derivative of the potential

`ss`: surface separation
`rab`: sum of nuclear radii
`b2y2a, b2y2b, b4y4a, b4y4b`: deformation effects through quadrupole and hexadecapole terms
`arg`: exponential term with diffuseness a0r
`ur`: nuclear potential with depth v0r

## <u>Function potent</u>

`rx`: radial distance
`v0`: total potential
`v1`: first derivative
`v2`: second derivative
`etak`: coulomb parameter (directly proportional to Z1Z2)

## <u>Function bar</u>

`rbar`: barrier radius
`vb`: barrier height
`homega`: curvature of barrier

`rmax`: maximum radius
`dr`: initial step size
`s`: direction indicator
`y1`: Temporary storage variable for finding the nuclear barrier height through a search method. It stores the previous potential value during an iterative process to locate the point where the potential function's sign changes.

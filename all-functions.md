# Formulas Used in Fortran Subroutines

Here are the formulas used in the Fortran subroutines `pot`, `potent`, and `bar`.

---

## 1. Subroutine `pot(rr, ur, dur, ddur)`

Computes the nuclear potential and its derivatives using a Woods-Saxon form factor with deformation corrections.

### Equations

1. **Deformation-adjusted radius**:  
   $$ ss = rr - rab - \underbrace{r_a \left(b2y2a + b4y4a\right)}_{\text{Projectile deformation}} - \underbrace{r_b \left(b2y2b + b4y4b\right)}_{\text{Target deformation}} $$  
   - $rab = r_a + r_b + 0.29$ (equilibrium separation)  
   - $b2y2a, b4y4a$: Quadrupole and hexadecapole deformation terms for the projectile  
   - $b2y2b, b4y4b$: Quadrupole and hexadecapole deformation terms for the target  

2. **Exponential argument**:  
   $$ \text{arg} = \exp\left(-\frac{ss}{a_{0r}}\right) $$  
   - $a_{0r}$: Diffuseness parameter (fixed at 0.63 fm)  

3. **Nuclear potential**:  
   $$ ur = -v_{0r} \cdot \frac{\text{arg}}{1 + \text{arg}} $$  
   - $v_{0r}$: Nuclear potential depth (adjusted from input `dv`)  
   $$v_{0r} = v_{0rx} * rred + dv - 20.0$$
   -  $v_{0rx}$: Base Nuclear Potential Term

$$v_{0rx} = 30.08 * (1. - 1.8 * (1. - 2.* \frac{wza}{wma}) * (1. - 2.* \frac{wzb}{wmb}))$$

   $wza, wzb$: Charges (atomic numbers) of the projectile and target, $wma, wmb$: Masses of the projectile and target.
   
   - $rred$: Reduced radius $$rred = \frac{r_{a}*r_{b}}{(r_{a}+r_{b})}$$where $r_a$ and $r_b$ are projectile/target radii.
   
   - $dv$: Input parameter (Typical range  `-10 < dv < 10`)

4. **First derivative of potential**:  
   $$ dur = -\frac{ur}{a_{0r} \cdot (1 + \text{arg})} $$  

5. **Second derivative of potential**:  
   $$ ddur = \frac{dur}{a_{0r}} \left(1 - \frac{2}{1 + \text{arg}}\right) $$  

---

## 2. Subroutine `potent(rx, v0, v1, v2)`

Computes the total potential (Coulomb + nuclear) and its first/second derivatives at radius $rx$.

### Equations

1. **Coulomb potential** (with deformation corrections):  
   $$ v0_{\text{coul}} = \frac{e_{\text{tak}}}{x} \left[1 + 0.6 \frac{b2y2a \cdot r_a^2}{x^2} + \frac{b4y4a \cdot r_a^4}{3x^4} + 0.6 \frac{b2y2b \cdot r_b^2}{x^2} + \frac{b4y4b \cdot r_b^4}{3x^4}\right] $$  
   - $e_{\text{tak}} = 1.43997 \cdot Z_p Z_t$ (Coulomb constant in MeV·fm)  
   - $x = rx$  

2. **First derivative of Coulomb potential**:  
   $$ v1_{\text{coul}} = -\frac{e_{\text{tak}}}{x^2} \left[1 + 1.8 \frac{b2y2a \cdot r_a^2}{x^2} + \frac{5b4y4a \cdot r_a^4}{3x^4} + 1.8 \frac{b2y2b \cdot r_b^2}{x^2} + \frac{5b4y4b \cdot r_b^4}{3x^4}\right] $$  

3. **Second derivative of Coulomb potential**:  
   $$ v2_{\text{coul}} = \frac{2e_{\text{tak}}}{x^3} \left[1 + 3.6 \frac{b2y2a \cdot r_a^2}{x^2} + \frac{15b4y4a \cdot r_a^4}{3x^4} + 3.6 \frac{b2y2b \cdot r_b^2}{x^2} + \frac{15b4y4b \cdot r_b^4}{3x^4}\right] $$  

4. **Total potential**:  
   $$ v0 = v0_{\text{coul}} + ur \quad (\text{from `pot' subroutine}) $$  
   $$ v1 = v1_{\text{coul}} + dur \quad (\text{first derivative}) $$  
   $$ v2 = v2_{\text{coul}} + ddur \quad (\text{second derivative}) $$  

---

## 3. Subroutine `bar(rbar, vb, homega)`

Finds the fusion barrier position ($r_{\text{bar}}$), height ($v_b$), and curvature ($\hbar\omega$).

### Algorithm

1. **Root-finding loop**:  
   - Iteratively adjust $r_{bz}$ until $v1$ (first derivative of potential) changes sign.  
   - Use bisection-like steps:       
   $$ \text{Initial: } r_{bz} = 20.0, \quad dr = 0.5 $$$$ \text{If } v1 \cdot y1 < 0 \implies \text{halve } dr \text{ and reverse search direction} $$  Terminate when $dr < 0.01$.  

2. **Barrier parameters**:  
   $$ r_{\text{bar}} = r_{bz} \quad (\text{barrier position}) $$  
   $$ v_b = v0 \quad (\text{barrier height}) $$  
   $$ \hbar\omega = hc \cdot \sqrt{\frac{-v2}{\text{redm} \cdot au}} \quad (\text{barrier curvature}) $$  
   - $hc = 197.3286$ MeV·fm ($\hbar c$)  
   - $\text{redm} = \frac{M_p M_t}{M_p + M_t}$: Reduced mass  
   - $au = 931.5016$ MeV/c² (atomic mass unit)  

---

## Key Variables

| Symbol           | Meaning                                  | Fortran Variable |
|------------------|------------------------------------------|------------------|
| $r_a, r_b$       | Projectile/Target radii                  | `ra`, `rb`       |
| $b2y2a, b4y4a$   | Projectile deformation terms             | Common block     |
| $e_{\text{tak}}$ | Coulomb constant                         | `etak`           |
| $a_{0r}$         | Nuclear potential diffuseness            | `a0r`            |
| $v_{0r}$         | Adjusted nuclear potential depth         | `v0r`            |

---

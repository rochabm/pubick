"""
Buckley-Leverett two-phase flow simulation
Oviedo, P. R. et al. (2026)

Schemes implemented:
  - FOU       : First-Order Upwind
  - QUICK     : Quadratic Upstream Interpolation for Convective Kinematics
  - SMART     : Sharp and Monotonic Algorithm for Realistic Transport
  - TOPUS     : Third-Order Polynomial Upwind Scheme
  - SOBUS     : (Bézier-based scheme)
  - HPUS      : High-order Polynomial Upwind Scheme
  - FDHPUS    : FD version of HPUS
  - ADBQUICKEST
  - CUBICK    : parameter set (0.5, 0.75)
  - CUBICK2   : parameter set (0.25, 0.45)
  - PUBICK    : parameter set (m1=3/10, m2=5/6)
  - PUBICK2   : parameter set (m1=0.493, m2=0.57)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from welge import exata_welge

import scienceplots
plt.style.use(['science','ieee'])

# =============================================================================
# Fractional flow helpers (Buckley-Leverett)
# =============================================================================

def f_bl(S, M):
    """Fractional flow function f(S) for mobility ratio c = M."""
    return S**2 / (S**2 + ((1 - S)**2)/M)

def df_bl(S, M):
    """Derivative of fractional flow: f'(S)"""
    return 2 * M * S * (1 - S) / (M * S**2 + (1 - S)**2)**2

# =============================================================================
# Initial condition
# =============================================================================

def initial_condition(x):
    u = np.zeros(len(x))
    u[x == 0] = 1.0
    return u

# =============================================================================
# Exact solution  (works for any mobility ratio c = M)
# =============================================================================

def exata(x, tf, c=1.0):
    """
    Exact Buckley-Leverett solution at time tf on grid x, f/ mobility ratio c=1.
    Uses Welge tangent construction for general c (welge.py).
    """
    f  = lambda S: f_bl(S, c)
    df = lambda S: df_bl(S, c)

    # Shock saturation: tangent from origin touches f(S)  =>  f'(Ss) = f(Ss)/Ss
    obj = lambda S: df(S) - f(S) / S
    Ss  = brentq(obj, 1e-6, 1.0 - 1e-6)
    dfs = df(Ss)   # shock speed

    N = len(x)
    u = np.zeros(N)
    for j in range(N):
        if x[j] <= 0:
            u[j] = 1.0
        else:
            xt = x[j] / tf
            if xt >= dfs:
                u[j] = 0.0
            else:
                # Rarefaction: invert df(S) = x/t
                try:
                    u[j] = brentq(lambda S: df(S) - xt, Ss, 1.0 - 1e-9)
                except ValueError:
                    u[j] = 1.0
    return u

# =============================================================================
# Shared flux-speed utility
# =============================================================================

def wave_speed(uL, uR, dt, dx, c):
    """Local wave speed (dt/dx) * df/du between cells uL and uR."""
    if uL == uR:
        return (dt / dx) * df_bl(uL, c)
    else:
        return (dt / dx) * (f_bl(uR, c) - f_bl(uL, c)) / (uR - uL)

def h(u, c):
    """Fractional flow evaluated at saturation u."""
    return f_bl(u, c)

# =============================================================================
# FOU – First-Order Upwind
# =============================================================================

def fou(u0, dt, dx, tf, c):
    nsteps = round(tf / dt)
    N = len(u0)
    u = u0.copy()
    unew = np.zeros(N)

    for _ in range(nsteps):
        for j in range(1, N - 1):
            v = wave_speed(u[j], u[j+1], dt, dx, c)
            top1 = h(u[j], c) if v > 0 else h(u[j+1], c)

            v = wave_speed(u[j-1], u[j], dt, dx, c)
            top2 = h(u[j-1], c) if v > 0 else h(u[j], c)

            unew[j] = u[j] - (dt / dx) * (top1 - top2)
        unew[0]  = u[0]
        unew[-1] = u[-1]
        u = unew.copy()
    return u

# =============================================================================
# Generic NVD-based upwind scheme
# A limiter function phi_func(phiU) -> normalised face value in [0,1]
# =============================================================================

def _nvd_scheme(u0, dt, dx, tf, c, limiter):
    """
    Generic NVD (Normalised Variable Diagram) finite-volume advection.
    limiter(phiU) returns the normalised face value phiF in [0,1].
    Falls back to upwind (phiF = phiU) outside [0,1].
    """
    nsteps = round(tf / dt)
    N = len(u0)
    u = u0.copy()
    unew = np.zeros(N)

    for _ in range(nsteps):
        for j in range(1, N - 1):

            # ---- right face (j+1/2) ----
            v = wave_speed(u[j], u[j+1], dt, dx, c)
            if v > 0:
                hU = h(u[j],   c)
                hD = h(u[j+1], c)
                hR = h(u[j-1], c)
            else:
                hU = h(u[j+1], c)
                hD = h(u[j],   c) if j <= N - 3 else hU
                hR = h(u[j+2], c) if j <= N - 3 else hU

            if hD == hR:
                top1 = hU
            else:
                phiU = (hU - hR) / (hD - hR)
                if 0.0 <= phiU <= 1.0:
                    top1 = hR + (hD - hR) * limiter(phiU)
                else:
                    top1 = hU

            # ---- left face (j-1/2) ----
            v = wave_speed(u[j-1], u[j], dt, dx, c)
            if v > 0:
                hU = h(u[j-1], c)
                hD = h(u[j],   c) if j >= 2 else hU
                hR = h(u[j-2], c) if j >= 2 else hU
            else:
                hU = h(u[j],   c)
                hD = h(u[j-1], c)
                hR = h(u[j+1], c)

            if hD == hR:
                top2 = hU
            else:
                phiU = (hU - hR) / (hD - hR)
                if 0.0 <= phiU <= 1.0:
                    top2 = hR + (hD - hR) * limiter(phiU)
                else:
                    top2 = hU

            unew[j] = u[j] - (dt / dx) * (top1 - top2)

        unew[0]  = u[0]
        unew[-1] = u[-1]
        u = unew.copy()
    return u

# =============================================================================
# QUICK
# =============================================================================

def quick(u0, dt, dx, tf, c):
    def lim(phi):
        if phi <= 5/6:
            return (6/8) * phi + 3/8
        else:
            return 1.0   # hD
    return _nvd_scheme(u0, dt, dx, tf, c, lim)

# =============================================================================
# SMART
# =============================================================================

def smart(u0, dt, dx, tf, c):
    def lim(phi):
        if phi < 1/6:
            return 3 * phi        # 3*hU - 2*hR in normalised coords = 3*phi
        elif phi <= 5/6:
            return (3 + 6*phi) / 8   # (3*hD + 6*hU - hR)/8 normalised
        else:
            return 1.0
    return _nvd_scheme(u0, dt, dx, tf, c, lim)

# =============================================================================
# TOPUS
# =============================================================================

def topus(u0, dt, dx, tf, c):
    def lim(phi):
        return 2*phi**4 - 3*phi**3 + 2*phi
    return _nvd_scheme(u0, dt, dx, tf, c, lim)

# =============================================================================
# SOBUS
# =============================================================================

def sobus(u0, dt, dx, tf, c):
    s3 = np.sqrt(3)
    def lim(phi):
        return (-(s3/2)*phi
                + (0.5 + s3/3) * (-(0.5*(3 - s3))
                + 3*np.sqrt((1/6)*(2 - s3) + (s3/3)*phi)))
    return _nvd_scheme(u0, dt, dx, tf, c, lim)

# =============================================================================
# HPUS
# =============================================================================

def hpus(u0, dt, dx, tf, c):
    def lim(phi):
        return phi * (-4*phi**4 + 10*phi**3 - 8*phi**2 + phi + 2)
    return _nvd_scheme(u0, dt, dx, tf, c, lim)

# =============================================================================
# FDHPUS  (a=1.5, b=0)
# =============================================================================

def fdhpus(u0, dt, dx, tf, c):
    a, b = 1.5, 0.0
    def lim(phi):
        return (4*(a+b-3)*phi**5 - 2*(6*a+4*b-17)*phi**4
                + (13*a+5*b-34)*phi**3 - (6*a+b-13)*phi**2 + a*phi)
    return _nvd_scheme(u0, dt, dx, tf, c, lim)

# =============================================================================
# ADBQUICKEST
# =============================================================================

def adbquickest(u0, dt, dx, tf, c, CFL):
    nsteps = round(tf / dt)
    N = len(u0)
    u = u0.copy()
    unew = np.zeros(N)

    numa = (1/3) - 0.5*abs(CFL) + (1/6)*CFL**2
    dena = 2 * ((7/12) - 0.5*CFL - 0.25*abs(CFL) + (1/6)*CFL**2)
    a    = numa / dena

    numb = (1/6)*CFL**2 - 2/3 - 0.5*abs(CFL) + CFL
    denb = 2 * ((-5/12) - 0.25*abs(CFL) + (1/6)*CFL**2 + 0.5*CFL)
    b    = numb / denb

    def lim(phi, A, B, C):
        if a <= phi <= b:
            return C + A + 0.5*(1 - abs(CFL))*(B - A) - (1/6)*(1 - CFL**2)*(B - 2*A)
        elif 0 <= phi < a:
            return C + (2 - CFL)*A
        elif b < phi <= 1:
            return C + B - B*CFL + A*CFL
        else:
            return C + A   # upwind (hU = C + A)

    for _ in range(nsteps):
        for j in range(1, N - 1):

            # right face
            v = wave_speed(u[j], u[j+1], dt, dx, c)
            if v > 0:
                hU_ = h(u[j],   c)
                hD_ = h(u[j+1], c)
                hR_ = h(u[j-1], c)
            else:
                if j <= N - 3:
                    hU_ = h(u[j+1], c)
                    hD_ = h(u[j],   c)
                    hR_ = h(u[j+2], c)
                else:
                    hU_ = h(u[j+1], c)
                    hD_ = hU_; hR_ = hU_

            if hD_ == hR_:
                top1 = hU_
            else:
                phiU = (hU_ - hR_) / (hD_ - hR_)
                A = hU_ - hR_; B = hD_ - hR_; C = hR_
                top1 = lim(phiU, A, B, C)

            # left face
            v = wave_speed(u[j-1], u[j], dt, dx, c)
            if v > 0:
                if j >= 2:
                    hU_ = h(u[j-1], c)
                    hD_ = h(u[j],   c)
                    hR_ = h(u[j-2], c)
                else:
                    hU_ = h(u[j-1], c); hD_ = hU_; hR_ = hU_
            else:
                hU_ = h(u[j],   c)
                hD_ = h(u[j-1], c)
                hR_ = h(u[j+1], c)

            if hD_ == hR_:
                top2 = hU_
            else:
                phiU = (hU_ - hR_) / (hD_ - hR_)
                A = hU_ - hR_; B = hD_ - hR_; C = hR_
                top2 = lim(phiU, A, B, C)

            unew[j] = u[j] - (dt / dx) * (top1 - top2)

        unew[0]  = u[0]
        unew[-1] = u[-1]
        u = unew.copy()
    return u

# =============================================================================
# CUBICK  (generic – accepts a1, b1, c1, d1)
# =============================================================================

def _cubick_core(u0, dt, dx, tf, c, a1, b1, c1, d1):
    A1 = (3*c1 - 6*a1) / (1 + 3*a1 - 3*c1)
    B1 = (3*a1)        / (1 + 3*a1 - 3*c1)
    C1 = (-1)          / (1 + 3*a1 - 3*c1)
    p  = (3*B1 - A1**2) / 3
    q1 = (2*A1**3 - 9*A1*B1) / 27

    def lim(phi):
        q  = q1 + C1*phi
        aq = np.sqrt((q**2)/4 + (p**3)/27)
        t  = (-q/2 + aq)**(1/3) - (q/2 + aq)**(1/3) - A1/3
        return 3*b1*t*(1-t)**2 + 3*d1*t**2*(1-t) + t**3

    return _nvd_scheme(u0, dt, dx, tf, c, lim)

def cubick(u0, dt, dx, tf, c):
    """CUBICK with parameters (a1=0.5, b1=0.75, c1=0.247622, d1=0.674287)."""
    return _cubick_core(u0, dt, dx, tf, c,
                        a1=0.5,  b1=0.75,
                        c1=0.247622138278052, d1=0.674286640691196)

def cubick2(u0, dt, dx, tf, c):
    """CUBICK with parameters (a1=0.25, b1=0.45, c1=0.256082, d1=0.735437)."""
    return _cubick_core(u0, dt, dx, tf, c,
                        a1=0.25, b1=0.45,
                        c1=0.256082165462314, d1=0.735437227905742)

# =============================================================================
# PUBICK  (generic – accepts m1, m2)
# =============================================================================

def _pubick_core(u0, dt, dx, tf, c, m1, m2):
    def lim(phi):
        if phi < 0.5:
            return ((3 / (4*(4*m1 - 1)**2))
                    * (4*m1*(4*m1 - 1)*phi
                       + (2*m1 - 1)*(2*m1 - np.sqrt(4*m1**2 + 2*(1 - 4*m1)*phi))))
        else:
            return (-(1 / (4*(4*m2 - 3)**2))
                    * ((4*m2 - 3)*(5 - 6*m2 + 4*(2 - 3*m2)*phi)
                       + (1 - 2*m2)*(1 - 2*m2
                                     + np.sqrt(4*m2**2 - 2 + 2*(3 - 4*m2)*phi))))

    return _nvd_scheme(u0, dt, dx, tf, c, lim)

def pubick(u0, dt, dx, tf, c):
    """PUBICK with parameters (m1=3/10, m2=5/6)."""
    return _pubick_core(u0, dt, dx, tf, c, m1=3/10, m2=5/6)

def pubick2(u0, dt, dx, tf, c):
    """PUBICK with parameters (m1=0.493, m2=0.57)."""
    return _pubick_core(u0, dt, dx, tf, c, m1=0.493, m2=0.57)

# =============================================================================
# Code to run BL problem and make plots
# =============================================================================

def run_buckley_leverett(M, tf=0.5, x0=0.0, xN=1.0, N=50, CFL=0.2,
                         filename=None, show=True):
    """
    Run all Buckley-Leverett schemes for a given mobility ratio M and plot results.

    Parameters
    ----------
    M        : float  Mobility ratio
    tf       : float  Final time
    x0, xN  : float  Domain bounds
    N        : int    Number of cells
    CFL      : float  CFL number
    filename : str    If given, save figure to this path (e.g. 'bl_M1.png')
    show     : bool   Whether to call plt.show()
    """
    dx = (xN - x0) / N
    dt = CFL * dx
    x  = np.linspace(x0, xN, N + 1)

    u0      = initial_condition(x)
    exact_w = exata_welge(x, tf, M)

    # Run all schemes
    solutions = {
        'CUBICK (0.25, 0.45)'                     : cubick2(u0, dt, dx, tf, M),
        r'PUBICK ($\mu_1$=3/10, $\mu_2$=5/6)'     : pubick (u0, dt, dx, tf, M),
        r'PUBICK ($\mu_1$=0.493, $\mu_2$=0.57)'   : pubick2(u0, dt, dx, tf, M),
        'TOPUS'                                   : topus  (u0, dt, dx, tf, M),
        'SMART'                                   : smart  (u0, dt, dx, tf, M),
        'QUICK'                                   : quick  (u0, dt, dx, tf, M),
    }

    styles = [
        ('royalblue',    '--'),
        ('crimson',      '-.'),
        ('seagreen',     (0, (5, 1))),
        ('darkorange',   ':'),
        ('mediumpurple', (0, (3, 1, 1, 1))),
        ('saddlebrown',  (0, (1, 1))),
    ]

    # Plot
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(x, exact_w, color='black', linestyle='-', linewidth=2.2, label='Exact (Welge)')

    for (label, sol), (color, ls) in zip(solutions.items(), styles):
        ax.plot(x, sol, color=color, linestyle=ls, linewidth=2, label=label)

    ax.set_xlabel('x', fontsize=14)
    ax.set_ylabel('u', fontsize=14)
    ax.set_ylim(-0.05, 1.05)
    ax.tick_params(labelsize=14)
    ax.legend(fontsize=12)
    # ax.set_title(f'Buckley-Leverett  (N={N}, CFL={CFL}, M={M}, tf={tf})', fontsize=13)
    ax.set_title(f'Buckley-Leverett (M={M})', fontsize=13)
    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=300)
        print(f"Saved: {filename}")
    if show:
        plt.show()

    return x, exact_w, solutions

# =============================================================================
# Main script 
# =============================================================================

if __name__ == "__main__":

    run_buckley_leverett(M=0.5, filename='buckley_leverett_M05.png')
    run_buckley_leverett(M=1.0, filename='buckley_leverett_M1.png')
    run_buckley_leverett(M=2.0, filename='buckley_leverett_M2.png')
    run_buckley_leverett(M=5.0, filename='buckley_leverett_M5.png')
    run_buckley_leverett(M=10.0, tf=0.25, filename='buckley_leverett_M10.png')

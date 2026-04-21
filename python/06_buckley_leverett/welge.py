"""
Exact Buckley-Leverett solution via Welge tangent construction.
Works for any mobility ratio M (= krw_end / kro_end).

Quick test: run this file directly.
    python exata_welge.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


# =============================================================================
# Core functions
# =============================================================================

def f_bl(S, M):
    """Fractional flow: f(S) = S^2 / (S^2 + (1-S)^2 / M)"""
    return S**2 / (S**2 + (1 - S)**2 / M)


def df_bl(S, M):
    """Derivative of fractional flow: f'(S)"""
    return 2 * M * S * (1 - S) / (M * S**2 + (1 - S)**2)**2


def welge_shock(M):
    """
    Find shock saturation Ss via Welge tangent construction.
    Solves:  f'(Ss) = f(Ss) / Ss   (tangent from origin touches f curve).
    Returns: Ss, shock speed (= f'(Ss))
    """
    obj = lambda S: df_bl(S, M) - f_bl(S, M) / S
    Ss  = brentq(obj, 1e-6, 1 - 1e-6)
    return Ss, df_bl(Ss, M)


def exata_welge(x, tf, M=1.0):
    """
    Exact Buckley-Leverett saturation profile at time tf.

    Parameters
    ----------
    x  : array-like   Spatial grid points
    tf : float        Final time
    M  : float        Mobility ratio (M = krw_end / kro_end). Default M=1.

    Returns
    -------
    u  : ndarray      Water saturation at each point in x
    """
    x = np.asarray(x, dtype=float)
    Ss, shock_speed = welge_shock(M)

    u = np.zeros_like(x)
    for j, xj in enumerate(x):
        if xj <= 0.0:
            u[j] = 1.0                          # inlet: fully saturated
        else:
            xt = xj / tf                        # characteristic speed x/t
            if xt >= shock_speed:
                u[j] = 0.0                      # ahead of shock
            else:
                # Rarefaction fan: invert f'(S) = x/t
                u[j] = brentq(lambda S: df_bl(S, M) - xt, Ss, 1 - 1e-9)

    return u


# =============================================================================
# Self-test
# =============================================================================

if __name__ == "__main__":

    x  = np.linspace(0, 1.5, 500)
    tf = 0.5

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # ---- Left: several mobility ratios ----
    ax = axes[0]
    M_values = [0.5, 1.0, 2.0, 5.0, 10.0]
    styles   = ['-', '--', '-.', ':', (0, (5, 1))]
    for M, ls in zip(M_values, styles):
        Ss, v_shock = welge_shock(M)
        u = exata_welge(x, tf, M)
        ax.plot(x, u, linestyle=ls, linewidth=2,
                label=f'M = {M}  (Ss={Ss:.3f}, vs={v_shock:.3f})')

    ax.set_xlabel('x', fontsize=13)
    ax.set_ylabel('S (water saturation)', fontsize=13)
    ax.set_title(f'Exact BL solution — various M  (tf={tf})', fontsize=13)
    ax.set_ylim(-0.05, 1.10)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # ---- Right: fractional flow curves + Welge tangent ----
    ax = axes[1]
    S_vec = np.linspace(0, 1, 300)
    colors = plt.cm.tab10(np.linspace(0, 0.5, len(M_values)))
    for M, col in zip(M_values, colors):
        f_vec = f_bl(S_vec, M)
        ax.plot(S_vec, f_vec, linewidth=1.8, color=col, label=f'M={M}')
        Ss, vs = welge_shock(M)
        fs = f_bl(Ss, M)
        # tangent line from origin through (Ss, fs)
        ax.plot([0, Ss], [0, fs], 'o--', color=col, markersize=5, linewidth=1)

    ax.set_xlabel('S', fontsize=13)
    ax.set_ylabel('f(S)', fontsize=13)
    ax.set_title('Fractional flow + Welge tangent (dot = shock point)', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('exata_welge.png', dpi=300)
    plt.show()

    # ---- Console report ----
    print(f"{'M':>6}  {'Ss':>8}  {'shock speed':>12}  {'f(Ss)':>8}")
    print("-" * 42)
    for M in M_values:
        Ss, vs = welge_shock(M)
        print(f"{M:>6.1f}  {Ss:>8.4f}  {vs:>12.4f}  {f_bl(Ss,M):>8.4f}")

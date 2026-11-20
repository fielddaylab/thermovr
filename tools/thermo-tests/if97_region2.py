"""Independent IAPWS-IF97 reference: Region 2 (superheated vapour) + Region 4 (saturation).

This module is an *oracle* for the ThermoVR isentropic-compression tests. It is
written straight from the IAPWS-IF97 release equations rather than ported from
Assets/Scripts/General/IAPWS97.cs, so a bug in the shipped backward functions
cannot be masked by the same bug in the oracle.

Native IF97 units are used throughout:
    p : MPa      T : K       v : m^3/kg
    h : kJ/kg    s : kJ/(kg K)

ThermoMath's public API uses Pa and J/(kg K); those conversions are asserted
explicitly by the test-suite, not hidden in here.
"""

import math

R = 0.461526      # kJ/(kg K), IF97 eq. 1
T_STAR = 540.0    # K, Region 2 reducing temperature
P_STAR = 1.0      # MPa, Region 2 reducing pressure

T_MIN = 273.15
T_MAX = 1073.15
P_MAX = 100.0
T_CRIT = 647.096
P_CRIT = 22.064

# ideal-gas part of the dimensionless Gibbs energy (IF97 Table 10): J_i, n_i
_IDEAL = (
    (0, -0.96927686500217e1),
    (1, 0.10086655968018e2),
    (-5, -0.56087911283020e-2),
    (-4, 0.71452738081455e-1),
    (-3, -0.40710498223928),
    (-2, 0.14240819171444e1),
    (-1, -0.43839511319450e1),
    (2, -0.28408632460772),
    (3, 0.21268463753307e-1),
)

# residual part of the dimensionless Gibbs energy (IF97 Table 11): I_i, J_i, n_i
_RESID = (
    (1, 0, -0.17731742473213e-2),
    (1, 1, -0.17834862292358e-1),
    (1, 2, -0.45996013696365e-1),
    (1, 3, -0.57581259083432e-1),
    (1, 6, -0.50325278727930e-1),
    (2, 1, -0.33032641670203e-4),
    (2, 2, -0.18948987516315e-3),
    (2, 4, -0.39392777243355e-2),
    (2, 7, -0.43797295650573e-1),
    (2, 36, -0.26674547914087e-4),
    (3, 0, 0.20481737692309e-7),
    (3, 1, 0.43870667284435e-6),
    (3, 3, -0.32277677238570e-4),
    (3, 6, -0.15033924542148e-2),
    (3, 35, -0.40668253562649e-1),
    (4, 1, -0.78847309559367e-9),
    (4, 2, 0.12790717852285e-7),
    (4, 3, 0.48225372718507e-6),
    (5, 7, 0.22922076337661e-5),
    (6, 3, -0.16714766451061e-10),
    (6, 16, -0.21171472321355e-2),
    (6, 35, -0.23895741934104e2),
    (7, 0, -0.59059564324270e-17),
    (7, 11, -0.12621808899101e-5),
    (7, 25, -0.38946842435739e-1),
    (8, 8, 0.11256211360459e-10),
    (8, 36, -0.82311340897998e1),
    (9, 13, 0.19809712802088e-7),
    (10, 4, 0.10406965210174e-18),
    (10, 10, -0.10234747095929e-12),
    (10, 14, -0.10018179379511e-8),
    (16, 29, -0.80882908646985e-10),
    (16, 50, 0.10693031879409),
    (18, 57, -0.33662250574171),
    (20, 20, 0.89185845355421e-24),
    (20, 35, 0.30629316876232e-12),
    (20, 48, -0.42002467698208e-5),
    (21, 21, -0.59056029685639e-25),
    (22, 53, 0.37826947613457e-5),
    (23, 39, -0.12768608934681e-14),
    (24, 76, 0.73087610595061e-28),
    (24, 81, 0.55414715350778e-16),
    (24, 84, -0.94369707241210e-16),
)


# --- Region 2 basic equation (dimensionless Gibbs energy and derivatives) ---

def _gamma_ideal(pi, tau):
    return math.log(pi) + sum(n * tau ** j for j, n in _IDEAL)


def _gamma_ideal_tau(pi, tau):
    return sum(n * j * tau ** (j - 1) for j, n in _IDEAL)


def _gamma_resid(pi, tau):
    return sum(n * pi ** i * (tau - 0.5) ** j for i, j, n in _RESID)


def _gamma_resid_pi(pi, tau):
    return sum(n * i * pi ** (i - 1) * (tau - 0.5) ** j for i, j, n in _RESID)


def _gamma_resid_tau(pi, tau):
    return sum(n * pi ** i * j * (tau - 0.5) ** (j - 1) for i, j, n in _RESID)


def v_pt(p, t):
    """Specific volume [m^3/kg] from pressure [MPa] and temperature [K]."""
    pi, tau = p / P_STAR, T_STAR / t
    # R is kJ/(kg K) so pressure must be expressed in kPa to yield m^3/kg
    return R * t / (p * 1000.0) * (1.0 + pi * _gamma_resid_pi(pi, tau))


def s_pt(p, t):
    """Specific entropy [kJ/(kg K)] from pressure [MPa] and temperature [K]."""
    pi, tau = p / P_STAR, T_STAR / t
    g = _gamma_ideal(pi, tau) + _gamma_resid(pi, tau)
    g_tau = _gamma_ideal_tau(pi, tau) + _gamma_resid_tau(pi, tau)
    return R * (tau * g_tau - g)


def h_pt(p, t):
    """Specific enthalpy [kJ/kg] from pressure [MPa] and temperature [K]."""
    pi, tau = p / P_STAR, T_STAR / t
    g_tau = _gamma_ideal_tau(pi, tau) + _gamma_resid_tau(pi, tau)
    return R * t * tau * g_tau


def u_pt(p, t):
    """Specific internal energy [kJ/kg] from pressure [MPa], temperature [K]."""
    return h_pt(p, t) - p * 1000.0 * v_pt(p, t)


# --- Region 4: saturation line (IF97 eq. 30 / 31) ---

_SAT = (
    0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
    0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
    -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
    0.65017534844798e3,
)


def psat_t(t):
    """Saturation pressure [MPa] for temperature [K]."""
    n = _SAT
    th = t + n[8] / (t - n[9])
    a = th * th + n[0] * th + n[1]
    b = n[2] * th * th + n[3] * th + n[4]
    c = n[5] * th * th + n[6] * th + n[7]
    return (2.0 * c / (-b + math.sqrt(b * b - 4.0 * a * c))) ** 4


def tsat_p(p):
    """Saturation temperature [K] for pressure [MPa]."""
    n = _SAT
    beta = p ** 0.25
    e = beta * beta + n[2] * beta + n[5]
    f = n[0] * beta * beta + n[3] * beta + n[6]
    g = n[1] * beta * beta + n[4] * beta + n[7]
    d = 2.0 * g / (-f - math.sqrt(f * f - 4.0 * e * g))
    return (n[9] + d - math.sqrt((n[9] + d) ** 2 - 4.0 * (n[8] + n[9] * d))) / 2.0


def is_superheated(p, t):
    """True when (p, T) lies in Region 2 (dry superheated vapour)."""
    if not (T_MIN <= t <= T_MAX) or not (0.0 < p <= P_MAX):
        return False
    if p >= P_CRIT:
        return t > T_CRIT
    return t > tsat_p(p)


# --- Backward functions solved from the forward equations ---
#
# At fixed pressure, Region 2 entropy is strictly increasing in temperature
# (ds/dT = cp/T > 0), so a bracketed bisection is both robust and exact to
# machine precision in a modest number of iterations.

def t_ps(p, s, tol=1e-12, max_iter=200):
    """Temperature [K] from pressure [MPa] and entropy [kJ/(kg K)] in Region 2."""
    lo = tsat_p(p) + 1e-9 if p < P_CRIT else T_CRIT
    hi = T_MAX
    s_lo, s_hi = s_pt(p, lo), s_pt(p, hi)
    if not (s_lo <= s <= s_hi):
        raise ValueError(
            "s=%r kJ/kgK out of Region 2 range [%r, %r] at p=%r MPa"
            % (s, s_lo, s_hi, p)
        )
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if s_pt(p, mid) < s:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol * max(1.0, mid):
            break
    return 0.5 * (lo + hi)


def v_ps(p, s):
    """Specific volume [m^3/kg] from pressure [MPa] and entropy [kJ/(kg K)]."""
    return v_pt(p, t_ps(p, s))


def isentropic_step(p1, t1, p2):
    """Walk an isentropic path from (p1, t1) to pressure p2.

    Returns (t2, v2, s) with s the conserved entropy in kJ/(kg K).
    """
    s = s_pt(p1, t1)
    t2 = t_ps(p2, s)
    return t2, v_pt(p2, t2), s

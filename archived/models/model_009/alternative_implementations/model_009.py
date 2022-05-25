from scipy.integrate import solve_ivp
from os.path import join, dirname, abspath
from mpmath import mp
import numpy as np

############################################################################################
# Model 009
#
# Alternative names: 'Full model (constant τS with infall and outflow)'
############################################################################################

# Constants

τS = 3.0  # [Gyr]
Twarm = 10_000.0  # [K]
Tcold = 100.0  # [K]
T1 = 50_000.0  # [K]
Σion = 8e-4  # [Mₒ pc^(-2)]
Σdiss = 1.5e-4  # [Mₒ pc^(-2)]
C2 = 0.074  # [Mₒ^2 pc^(-4) Gyr]
C4 = 798e-3  # [Mₒ^2 pc^(-4) Gyr]
T = 13.7  # [Gyr]
ΣI = 100.0  # [Mₒ pc^(-2)]
τ = 2.0  # [Gyr] - {10^(-3), 2, 4, 10^3}
α = 0.3
ηi_lim = 1000.0
ηd_lim = 75.0
ω = 2.0
Zsun = 0.006
Zeff = 1e-3 * Zsun
Zsn = 0.06
R = 0.17

# System of equations

"""
	Ionized gas density:       i(t) -> y[0]
	Atomic gas density:        a(t) -> y[1]
	Molecular gas density:     m(t) -> y[2]
	Metal density:             z(t) -> y[3]
	Stellar density:           s(t) -> y[4]		

    Each equation has Mₒ pc^(-2) Gyr^(-1) [solar_mass * parsec^(-2) * years^(-9)]
	as its units on the LHS and RHS
"""


def system(t, y):

    # Variables
    i, a, m, z, s = y

    # Auxiliary equations
    g = i + a + m
    tot = g + s
    star_elem = m + α * a
    ψ = star_elem / τS
    Ii = ΣI * mp.exp(-t / τ) / (τ * (1.0 - mp.exp(-T / τ)))
    outflow = ω * ψ / g
    Oi = outflow * i
    Oa = outflow * a
    Om = outflow * m
    Oz = outflow * z
    ηion = ηi_lim * (1 - mp.exp(-a / Σion))
    ηdiss = ηd_lim * (1 - mp.exp(-m / Σdiss))
    Z = z / g
    τR = C2 * (1 + T1 * ψ / (Twarm * g)) / (g * tot)
    τC = C4 * (1 + T1 * ψ / (Tcold * g)) * Zsun / (g * tot * (Z + Zeff))
    recombination = i / τR
    cloud_formation = a / τC

    # ODE system
    return [
        Ii - Oi - recombination + (ηion + R) * ψ,
        -Oa - cloud_formation + recombination + (ηdiss - ηion - α * a / star_elem) * ψ,
        -Om + cloud_formation - (ηdiss + m / star_elem) * ψ,
        -Oz + (Zsn * R - Z) * ψ,
        (1 - R) * ψ,
    ]


############################################################################################
# Solving the system
############################################################################################

# Integration span

t_start = 0.0  # [Gyr]
t_end = 1.0  # [Gyr]
t_span = [t_start, t_end]

# Initial conditions

g0 = 1e-3 * 14.8285  # [Mₒ pc^(-2)]
i0 = g0 * 0.6  # [Mₒ pc^(-2)]
a0 = g0 * 0.2  # [Mₒ pc^(-2)]
m0 = g0 * 0.2  # [Mₒ pc^(-2)]
z0 = g0 * 1e-4  # [Mₒ pc^(-2)]
s0 = g0 * 0.0  # [Mₒ pc^(-2)]
init_cond = [i0, a0, m0, z0, s0]

# Parameters

# There are no parameters in this model.

params = None

# Integration

sol_BDF = solve_ivp(
    fun=system,
    t_span=t_span,
    y0=init_cond,
    method="BDF",
    dense_output=True,
    rtol=1e-6,
    atol=1e-8,
    first_step=(t_end - t_start) * 1e-5,
    args=params,
)
t = np.logspace(start=-5, stop=np.log10(t_end), num=int(1e4))
y = sol_BDF.sol(t)

# Save data in file

data = np.stack((t, y[0], y[1], y[2], y[3], y[4]), axis=1)
np.savetxt(
    join(dirname(abspath(__file__)), "plots/python.dat"),
    data,
    header="t \t i(t) \t a(t) \t m(t) \t z(t) \t s(t) \n",
)

num = y[4][-1]
ndigits = 3 - (int(mp.floor(mp.log10(abs(num)))) + 1)  # Only works for numbers < 1.
print("Stellar density after", t_end, "Gyr:", round(num, ndigits), "Mₒpc^(-2)")

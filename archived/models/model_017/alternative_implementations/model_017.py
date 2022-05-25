from scipy.integrate import solve_ivp
from os.path import join, dirname, abspath
from mpmath import mp
import numpy as np

############################################################################################
# Model 017
#
# Alternative names: 'New full model', 'Model D'
############################################################################################

# Constants

K = 1 / 400.0  # [Mₒ^(-2/3) pc^(4/3) Gyr^(-1)]
C1 = 3.03e-6  # [Gyr]
C2 = 8.74e-5  # [Gyr]
ηion = 955.29
ηdiss = 380.93
Zsun = 0.0134
Zeff = 1e-3 * Zsun
Zsn = 0.09
R = 0.18

# System of equations

"""
	Ionized gas density:       i(t) -> y[0]
	Atomic gas density:        a(t) -> y[1]
	Molecular gas density:     m(t) -> y[2]
	Metal density:             z(t) -> y[3]
	Stellar density:           s(t) -> y[4]		

    Each equation has Mₒ pc^(-3) Gyr^(-1) [solar_mass * parsec^(-3) * years^(-9)]
	as its units on the LHS and RHS
"""


def system(t, y):

    # Variables
    i, a, m, z, s = y

    # Auxiliary equations
    g = i + a + m
    tot = g + s
    τS = m ** (1 / 3.0) / (K * g ** (1 / 3.0) * tot ** (2 / 3.0))
    ψ = m / τS
    Z = z / g
    τR = C1 / i
    τC = C2 / ((a + m) * (Z + Zeff))
    recombination = i / τR
    cloud_formation = a / τC

    # ODE system
    return [
        -recombination + (ηion + R) * ψ,
        -cloud_formation + recombination + (ηdiss - ηion) * ψ,
        cloud_formation - (ηdiss + 1) * ψ,
        (Zsn * R - Z) * ψ,
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

g0 = 1e-3 * 14.8285  # [Mₒ pc^(-3)]
i0 = g0 * 0.6  # [Mₒ pc^(-3)]
a0 = g0 * 0.2  # [Mₒ pc^(-3)]
m0 = g0 * 0.2  # [Mₒ pc^(-3)]
z0 = g0 * 1e-4  # [Mₒ pc^(-3)]
s0 = g0 * 0.0  # [Mₒ pc^(-3)]
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

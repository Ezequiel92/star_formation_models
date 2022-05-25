from scipy.integrate import solve_ivp
from os.path import join, dirname, abspath
from mpmath import mp
import numpy as np

############################################################################################
# Model 014
#
# Alternative names: 'Simplified model'
############################################################################################

# Constants

K = 1 / 400.0  # [Mₒ^(-2/3) pc^(4/3) Gyr^(-1)]
Twarm = 10_000.0  # [K]
Tcold = 100.0  # [K]
T1 = 50_000.0  # [K]
Σion = 8e-4  # [Mₒ pc^(-2)]
Σdiss = 1.5e-4  # [Mₒ pc^(-2)]
C2 = 0.074  # [Mₒ^2 pc^(-4) Gyr]
C4 = 798e-3  # [Mₒ^2 pc^(-4) Gyr]
α = 7.0
ηi_lim = 955.29
ηd_lim = 380.93
Zsun = 0.0134
Zeff = 1e-3 * Zsun
Zsn = 0.2
R = 0.17

# System of equations

"""
    Ionized gas fraction:       i(t) / g -> y[0]
    Atomic gas fraction:        a(t) / g -> y[1]
    Molecular gas fraction:     m(t) / g -> y[2]
    Metal fraction:             z(t) / g -> y[3]
    Stellar fraction:           s(t) / g -> y[4]		
     
    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
"""


def system(t, y, g):

    # Variables
    i, a, m, z, s = y

    # Auxiliary equations
    star_elem = m + α * a
    ψ = K * (star_elem * g) ** (2 / 3)
    ηion = ηi_lim * (1 - mp.exp(-g * a / Σion))
    ηdiss = ηd_lim * (1 - mp.exp(-g * m / Σdiss))
    Z = z
    τR = C2 * (1 + T1 * ψ / Twarm) / (g * g)
    τC = C4 * (1 + T1 * ψ / Tcold) * Zsun / (g * g * (Z + Zeff))
    recombination = i / τR
    cloud_formation = a / τC

    # ODE system
    return [
        -recombination + (ηion + R) * ψ,
        -cloud_formation + recombination + (ηdiss - ηion - α * a / star_elem) * ψ,
        cloud_formation - (ηdiss + m / star_elem) * ψ,
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

i0 = 0.6
a0 = 0.2
m0 = 0.2
z0 = 1e-4
s0 = 0.0
init_cond = [i0, a0, m0, z0, s0]

# Parameters

params = (1.0,)  # g [Mₒ pc^(-2)]

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
print("Stellar fraction after", t_end, "Gyr:", round(num, ndigits))

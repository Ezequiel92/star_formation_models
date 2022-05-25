from scipy.integrate import solve_ivp
import numpy as np
import math as mt
from os.path import join, dirname, abspath

############################################################################################
# Model 019
############################################################################################

# Constants

C1 = 0.809  # [Gyr Mₒ^(1/2) pc^(-3/2)]
C2 = 3.024e-6  # [Gyr Mₒ pc^(-3)]
C3 = 8.745e-5  # [Gyr Mₒ pc^(-3)]
Zsun = 0.0134
Zeff = 1e-3 * Zsun
Zsn = 0.09
ηion = 955.29
ηdiss = 380.93
R = 0.18

# System of equations

"""
    Ionized gas fraction:       i(t) / ρ -> y[1]
    Atomic gas fraction:        a(t) / ρ -> y[2]
    Molecular gas fraction:     m(t) / ρ -> y[3]
    Metal fraction:             z(t) / ρ -> y[4]
    Stellar fraction:           s(t) / ρ -> y[5]	

    where ρ = i(t) + a(t) + m(t) + s(t)
     
    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
"""


def system(t, y, ρ0, g0):

    # Variables
    i, a, m, z, s = y

    # Auxiliary equations
    g = i + a + m
    τS = (C1 * g0) / mt.sqrt(g * ρ0)
    τR = C2 / (i * ρ0)
    τC = C3 / (g * ρ0 * (z + Zeff))
    recombination = i / τR
    cloud_formation = a / τC
    ψ = m / τS

    # ODE system
    return [
        -recombination + (ηion + R) * ψ,
        -cloud_formation + recombination + (ηdiss - ηion) * ψ,
        cloud_formation - (1 + ηdiss) * ψ,
        (Zsn * R - z) * ψ,
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

i0 = 0.9
a0 = 0.05
m0 = 0.05
z0 = 1e-4
s0 = 0.0
init_cond = [i0, a0, m0, z0, s0]

# Parameters

params = (0.9 / 40.46, i0 + a0 + m0)

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
ndigits = 3 - (int(mt.floor(mt.log10(abs(num)))) + 1)  # Only works for numbers < 1.
print("Stellar fraction after", t_end, "Gyr:", round(num, ndigits))

from scipy.integrate import solve_ivp
import numpy as np
import math as mt
from os.path import join, dirname, abspath

############################################################################################
# Model 018
#
# Alternative names: 'Fully coupled model (constant τS)', 'Model E'
############################################################################################

# Constants

τS = 2.6  # [Gyr]
C1 = 3.03e-6  # [Gyr]
C2 = 8.74e-5  # [Gyr]
P = 40.4  # [cm^(-3) Mₒ^(-1) pc^3]
Zsun = 0.0134
Zeff = 1e-3 * Zsun
Zsn = 0.09
ηion = 955.29
ηdiss = 380.93
R = 0.18

# System of equations

"""
	Ionized gas fraction:       i(t) / tot -> y[1]
    Atomic gas fraction:        a(t) / tot -> y[2]
    Molecular gas fraction:     m(t) / tot -> y[3]
    Metal fraction:             z(t) / tot -> y[4]
    Stellar fraction:           s(t) / tot -> y[5]	

    where tot = i(t) + a(t) + m(t) + z(t) + s(t)
     
    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
"""


def system(t, y, tot0):

    # Variables
    i, a, m, z, s = y

    # Auxiliary equations
    g = i + m + a + z
    ψ = m / τS
    ig = i / g
    ag = a / g
    mg = m / g
    zg = z / g
    τR = C1 / (i * tot0)
    τC = C2 / ((a + m) * tot0 * (zg + Zeff))
    recombination = i / τR
    cloud_formation = a / τC

    # ODE system
    return [
        -recombination + (ηion + (1 - Zsn) * R - ig) * ψ,
        -cloud_formation + recombination + (ηdiss - ηion - ag) * ψ,
        cloud_formation - (ηdiss + mg) * ψ,
        (Zsn * R - zg) * ψ,
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

params = (1e-3 / P,)  # tot0 [Mₒ pc^(-3)] = n / P

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

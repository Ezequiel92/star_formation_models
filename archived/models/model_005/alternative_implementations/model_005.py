from scipy.integrate import solve_ivp
import numpy as np
import math as mt
from os.path import join, dirname, abspath

############################################################################################
# Model 005
#
# Alternative names: 'Modified basic model with a non-constant number density'
############################################################################################

# Constants

K1 = 0.001656  # [pc^4 Mₒ^(-2) cm^(-3)]
K2 = 0.00984  # [pc^4 Mₒ^(-2) cm^(-3)]
τS = 2.6  # [Gyr]
Rh = 1.9  # [cm^3 Gyr^-1]
σνB = 8158.0  # [cm^3 Gyr^-1]
Zsun = 0.0134
Zeff = 1e-3 * Zsun
Zsn = 0.09
ηion = 955.29
ηdiss = 380.93
R = 0.18

# System of equations

"""
	Ionized gas fraction:       i(t) / g -> y[0]
    Atomic gas fraction:        a(t) / g -> y[1]
    Molecular gas fraction:     m(t) / g -> y[2]
    Metal fraction:             z(t) / g -> y[3]
    Stellar fraction:           s(t) / g -> y[4]		

    where g = i(t) + a(t) + m(t)
     
    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
"""


def system(t, y, g, Z):

    # Variables
    i, a, m, z, s = y

    # Auxiliary equations
    ne = K1 * (1 + z) * g * g
    nh = K2 * (1 + z) * g * g
    τR = 1 / (ne * σνB)
    τC = Zsun / (2 * nh * Rh * (Z + Zeff))
    recombination = i / τR
    cloud_formation = a / τC
    ψ = m / τS

    # ODE system
    return [
        -recombination + ηion * ψ,
        -cloud_formation + recombination + (ηdiss - ηion) * ψ,
        cloud_formation - (1 + ηdiss) * ψ,
        (Zsn * R - Z) * ψ,
        ψ,
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

params = (
    1e-3 * 14.8285,  # g [Mₒ pc^(-2)]
    z0,  # Z
)

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

using DifferentialEquations
using DelimitedFiles

############################################################################################
# Model 019
############################################################################################

# Constants

c₁ = 0.809       # [Gyr Mₒ^(1/2) pc^(-3/2)]
c₂ = 3.024e-6    # [Gyr Mₒ pc^(-3)]
c₃ = 8.745e-5    # [Gyr Mₒ pc^(-3)]
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
function system!(dydt, y, p, t)

    # Variables
    i, a, m, z, s = y

    # Parameters
    ρ₀, g₀ = p

    # Auxiliary equations
    g = i + a + m
    τS = (c₁ * g₀) / sqrt(g * ρ₀)
    τR = c₂ / (i * ρ₀)
    τC = c₃ / (g * ρ₀ * (z + Zeff))
    recombination = i / τR
    cloud_formation = a / τC
    ψ = m / τS

    # ODE system
    dydt[1] = -recombination + (ηion + R) * ψ
    dydt[2] = -cloud_formation + recombination + (ηdiss - ηion) * ψ
    dydt[3] = cloud_formation - (1 + ηdiss) * ψ
    dydt[4] = (Zsn * R - z) * ψ
    dydt[5] = (1 - R) * ψ

end

############################################################################################
# Solving the system
############################################################################################

# Integration span

t_start = 0.0    # [Gyr]
t_end = 1.0      # [Gyr]
tspan = (t_start, t_end)

# Initial conditions

i₀ = 0.9
a₀ = 0.05
m₀ = 0.05
z₀ = 1e-4
s₀ = 0.0
ic = [i₀, a₀, m₀, z₀, s₀]

# Parameters

params = [
    0.9 / 40.46,     # ρ₀ [Mₒ^(-1) pc^(3)] = n [cm^(-3)] / 40.3
    i₀ + a₀ + m₀,    # g₀
]

# Integration

lf = log10(t_end)
li = lf - 5.0
points = 1e4
step = 5.0 / points
times = [10^i for i = li:step:lf]

prob = ODEProblem(system!, ic, tspan, params)
sol = solve(prob, reltol = 1e-6, abstol = 1e-8, saveat = times)

# Save data in file

open(joinpath(@__DIR__, "plots/julia.dat"), "w") do file
    write(file, "t \t i(t) \t a(t) \t m(t) \t z(t) \t s(t) \n")
    writedlm(file, zip(sol.t, sol[1, :], sol[2, :], sol[3, :], sol[4, :], sol[5, :]))
end

println(
    "Stellar fraction after ",
    t_end,
    "Gyr: ",
    round(sol(t_end, idxs = 5), sigdigits = 3),
)

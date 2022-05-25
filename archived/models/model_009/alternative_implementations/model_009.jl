using DifferentialEquations
using DelimitedFiles

############################################################################################
# Model 009
#
# Alternative names: 'Full model (constant τS with infall and outflow)'
############################################################################################

# Constants

τS = 3.0            # [Gyr]
Twarm = 10_000.0    # [K]
Tcold = 100.0       # [K]
T1 = 50_000.0       # [K]
Σion = 8e-4         # [Mₒ pc^(-2)]
Σdiss = 1.5e-4      # [Mₒ pc^(-2)]
C2 = 0.074          # [Mₒ^2 pc^(-4) Gyr]
C4 = 798e-3         # [Mₒ^2 pc^(-4) Gyr]
T = 13.7            # [Gyr]
ΣI = 100.0          # [Mₒ pc^(-2)]
τ = 2.0         	# [Gyr] - {10^(-3), 2, 4, 10^3}
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
	Ionized gas density:       i(t) -> y[1]
	Atomic gas density:        a(t) -> y[2]
	Molecular gas density:     m(t) -> y[3]
	Metal density:             z(t) -> y[4]
	Stellar density:           s(t) -> y[5]		
     
    Each equation has Mₒ pc^(-2) Gyr^(-1) [solar_mass * parsec^(-2) * years^(-9)]
	as its units on the LHS and RHS
"""
function system!(dydt, y, p, t)

    # Variables
	i, a, m, z, s = y
	
	# Auxiliary equations
    g = i + a + m
    tot = g + s
	star_elem = m + α * a
	ψ = star_elem / τS
    Ii = ΣI * exp(-t / τ) / (τ * (1.0 - exp(-T / τ)))
    outflow = ω * ψ / g
    Oi = outflow * i
    Oa = outflow * a
    Om = outflow * m
    Oz = outflow * z
    ηion = ηi_lim * (1 - exp(-a / Σion))
    ηdiss = ηd_lim * (1 - exp(-m / Σdiss))
    Z = z / g
    τR = C2 * (1 + T1 * ψ / (Twarm * g)) / (g * tot)
    τC = C4 * (1 + T1 * ψ / (Tcold * g)) * Zsun / (g * tot * (Z + Zeff))
    recombination = i / τR
    cloud_formation = a / τC
	
	# ODE system
	dydt[1] = Ii - Oi - recombination + (ηion + R) * ψ
    dydt[2] = -Oa - cloud_formation + recombination + (ηdiss - ηion - α * a / star_elem) * ψ
    dydt[3] = -Om + cloud_formation - (ηdiss + m / star_elem) * ψ
	dydt[4] = -Oz + (Zsn * R - Z) * ψ
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

g₀ = 1e-3 * 14.8285    # [Mₒ pc^(-2)]
i₀ = g₀ * 0.6          # [Mₒ pc^(-2)]
a₀ = g₀ * 0.2          # [Mₒ pc^(-2)]
m₀ = g₀ * 0.2          # [Mₒ pc^(-2)]
z₀ = g₀ * 1e-4         # [Mₒ pc^(-2)] 
s₀ = g₀ * 0.0          # [Mₒ pc^(-2)]
init_cond = [i₀, a₀, m₀, z₀, s₀]

# Parameters

# There are no parameters in this model.

params = nothing

# Integration

log_start = -5
log_end = log10(t_end)
points = 1e4
times = [10^i for i in log_start:((log_end - log_start) / points):log_end]

prob = ODEProblem(system!, init_cond, tspan)
sol = solve(prob, Rodas5(), reltol = 1e-6, abstol = 1e-8, saveat = times)

# Save data in file

open(joinpath(@__DIR__, "plots/julia.dat"), "w") do file
    write(file, "t \t i(t) \t a(t) \t m(t) \t z(t) \t s(t) \n")
    writedlm(file, zip(sol.t, sol[1,:], sol[2,:], sol[3,:], sol[4,:], sol[5,:]))
end

println(
    "Stellar density after ", 
    t_end, 
    "Gyr: ", 
    round(sol(t_end, idxs = 5), sigdigits = 3), 
    " Mₒpc^(-2)",
)

using DifferentialEquations
using DelimitedFiles

############################################################################################
# Model 007
#
# Alternative names: 'Model B', 'Ascasibar et al. model (constant τS)'
############################################################################################

# Constants

C1 = 0.074      		# [Mₒ^2 pc^(-4) Gyr]
C2 = 360e-3     		# [Mₒ^2 pc^(-4) Gyr]
τS = 2.6                # [Gyr]
Zsun = 0.0134
Zeff = 1e-3 * Zsun
Zsn = 0.09
ηion = 955.29
ηdiss = 380.93
R = 0.18

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
    Z = z / g
    τR = C1 / (g * tot)
    τC = C2 / (g * tot * (Z + Zeff))
    recombination = i / τR
    cloud_formation = a / τC
	ψ = m / τS
	
	# ODE system
    dydt[1] = -recombination + (ηion + R) * ψ         
    dydt[2] = -cloud_formation + recombination + (ηdiss - ηion) * ψ
    dydt[3] = cloud_formation - (1 + ηdiss) * ψ
	dydt[4] = (Zsn * R - Z) * ψ
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

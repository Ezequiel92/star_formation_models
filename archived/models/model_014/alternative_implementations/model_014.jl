using DifferentialEquations
using DelimitedFiles

############################################################################################
# Model 014
#
# Alternative names: 'Simplified model'
############################################################################################

# Constants

K = 1 / 400.0       # [Mₒ^(-2/3) pc^(4/3) Gyr^(-1)]
Twarm = 10_000.0    # [K]
Tcold = 100.0       # [K]
T1 = 50_000.0       # [K]
Σion = 8e-4         # [Mₒ pc^(-2)]
Σdiss = 1.5e-4      # [Mₒ pc^(-2)]
C2 = 0.074          # [Mₒ^2 pc^(-4) Gyr]
C4 = 798e-3         # [Mₒ^2 pc^(-4) Gyr]
α = 7.0
ηi_lim = 955.29
ηd_lim = 380.93
Zsun = 0.0134
Zeff = 1e-3 * Zsun
Zsn = 0.2
R = 0.17

# System of equations

"""
    Ionized gas fraction:       i(t) / g -> y[1]
    Atomic gas fraction:        a(t) / g -> y[2]
    Molecular gas fraction:     m(t) / g -> y[3]
    Metal fraction:             z(t) / g -> y[4]
    Stellar fraction:           s(t) / g -> y[5]		
     
    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
"""
function system!(dydt, y, p, t)

    # Variables
	i, a, m, z, s = y

    # Parameters
	g = params[1]
	
	# Auxiliary equations
	star_elem = m + α * a
    ψ = K * (star_elem * g)^(2 / 3)
    ηion = ηi_lim * (1 - exp(-g * a / Σion))
    ηdiss = ηd_lim * (1 - exp(-g * m / Σdiss))
	Z = z
    τR = C2 * (1 + T1 * ψ / Twarm) / (g * g)
    τC = C4 * (1 + T1 * ψ / Tcold) * Zsun / (g * g * (Z + Zeff))
    recombination = i / τR
    cloud_formation = a / τC
	
	# ODE system
	dydt[1] = -recombination + (ηion + R) * ψ
    dydt[2] = -cloud_formation + recombination + (ηdiss - ηion - α * a / star_elem) * ψ
    dydt[3] = cloud_formation - (ηdiss + m / star_elem) * ψ
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

i₀ = 0.6    
a₀ = 0.2   
m₀ = 0.2  
z₀ = 1e-4  
s₀ = 0.0
init_cond = [i₀, a₀, m₀, z₀, s₀]

# Parameters

params = [
    1.0,    # g [Mₒ pc^(-2)]
]     

# Integration

log_start = -5
log_end = log10(t_end)
points = 1e4
times = [10^i for i in log_start:((log_end - log_start) / points):log_end]

prob = ODEProblem(system!, init_cond, tspan)
sol = solve(prob, ABDF2(), reltol = 1e-8, abstol = 1e-8, saveat = times)

# Save data in file

open(joinpath(@__DIR__, "plots/julia.dat"), "w") do file
    write(file, "t \t i(t) \t a(t) \t m(t) \t z(t) \t s(t) \n")
    writedlm(file, zip(sol.t, sol[1,:], sol[2,:], sol[3,:], sol[4,:], sol[5,:]))
end

println(
    "Stellar fraction after ", 
    t_end, 
    "Gyr: ", 
    round(sol(t_end, idxs = 5), sigdigits = 3),
)

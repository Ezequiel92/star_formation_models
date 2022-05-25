using DifferentialEquations
using DelimitedFiles

############################################################################################
# Model 004
#
# Alternative names: 'Modified basic model with mass recycling'
############################################################################################

# Constants

τS = 2.6        # [Gyr]
Rh = 1.9        # [cm^3 Gyr^-1]
σνB = 8158.0    # [cm^3 Gyr^-1]
Zsun = 0.0134
Zeff = 1e-3 * Zsun
Zsn = 0.09
ηion = 955.29
ηdiss = 380.93
R = 0.18

# System of equations

"""
	Ionized gas fraction:       i(t) / g -> y[1]
    Atomic gas fraction:        a(t) / g -> y[2]
    Molecular gas fraction:     m(t) / g -> y[3]
    Metal fraction:             z(t) / g -> y[4]
    Stellar fraction:           s(t) / g -> y[5]	

    where g = i(t) + a(t) + m(t)
     
    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
"""
function system!(dydt, y, p, t)

    # Variables
	i, a, m, z, s = y
	
	# Parameters
	n, Z = params
	
	# Auxiliary equations
    τR = 1 / (n * σνB)
    τC = Zsun / (2 * n * Rh * (Z + Zeff))
    recombination = i / τR
    cloud_formation = a / τC
    ψ = m / τS
	
	# ODE system
    dydt[1] = -recombination + ηion * ψ + R * ψ          
    dydt[2] = -cloud_formation + recombination + (ηdiss - ηion) * ψ
    dydt[3] = cloud_formation - (1 + ηdiss) * ψ
    dydt[4] = (Zsn * R - Z) * ψ
	dydt[5] = ψ - R * ψ 
	
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
	1e-3,    # n [cm^(-3)]
	z₀,      # Z
]  

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
    "Stellar fraction after ", 
    t_end, 
    "Gyr: ", 
    round(sol(t_end, idxs = 5), sigdigits = 3),
)

### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 8ae23d20-81ad-11eb-0644-4796c8c1e8c0
using Unitful, UnitfulAstro, UnitfulMoles, PlutoUI

# ╔═╡ b8930ef2-e0a3-4ff9-9cce-88dcaf486cef
TableOfContents(title="🔢 Constants")

# ╔═╡ bf6ad952-81b0-11eb-1de7-57344a180233
md"### Recombination rate, Carlos $[\text{cm}^3 \,\text{s}^{-1}]$ $\rightarrow$ me $[\text{cm}^3 \, \text{Gyr}^{-1}]$"

# ╔═╡ 44bcd24e-81ae-11eb-392f-43c49f59f6da
f(T) = (4.1e-10 * T^(-0.8))u"cm^3/s";

# ╔═╡ 21d3bdb0-81c4-11eb-3e58-973170b3ce20
4.1e-10u"cm^3/s" |> u"cm^3" * UnitfulAstro.Gyr^-1

# ╔═╡ f68db520-81c3-11eb-0463-412bda3df6c6
md"### Constants and conversion factors"

# ╔═╡ 47281080-81af-11eb-0bc3-8762ecd9a5f0
begin
	
	const Tc = 100u"K"	                        # T_cold
	const Tw = 1e4u"K"	                        # T_warm
	const Teff = 1e4u"K"	                    # T_eff
	const Zsun = 0.0134                         # Solar metallicity
	const σνB = f(1e4)
	const galaxythickness = 0.6UnitfulAstro.kpc
	
	# From `out` file of run_A_02
	# It should be the same in every other isolated simulation
	const a_scale = 1.0
	const a3inv = 1/ a_scale^3
	
	const SEC_PER_YEAR = 3.15576e7              # 1.0 yr in cm
	const SOLAR_MASS = 1.989e33                 # 1.0 M⊙ in g
	const PARSEC = 3.085678e18                  # 1.0 pc in cm
	const PROTONMASS = 1.67262178e-24           # 1.0 mp in g
	const HYDROGEN_MASSFRAC = 0.76              # Dimensionless
	const HUBBLE = 3.2407789e-18				# 1.0 H (Hubble parameter) in h/sec
	const GRAVITY = 6.672e-8                    # 1.0 G (gravitational constant) in cm^3 g^-1 s^-2
	
	# Cosmological constants in inernal units
	const Hubble_internal_units = 0.1
	const G_internal_units = 43007.1
	
	# Cosmological dimensionless parameters
	const Omega0 = 0.0
	const OmegaLambda = 0.0
	const OmegaBaryon = 0.0
	const HubbleParam = 1.0
	
	# Conversion factors from internal units to cgs
	const UnitMass_in_g = 1.989e+43             # 1.0e10 solar masses in g
	const UnitTime_in_s = 3.08568e+16           # 1.0 s*kpc/km in s
	const UnitVelocity_in_cm_per_s = 100000     # 1 km/s in cm/s
	const UnitDensity_in_cgs = 6.76991e-22      # 1.0e10 M⊙/kpc^3 in g/cm^3
	const UnitEnergy_in_cgs = 1.989e+53         # 1.0e10 M⊙*km^2/s^2 in g*cm^2/s^2 (ergios)
	const UnitLength_in_cm = 3.08568e+21        # 1.0 kpc in cm
	
	# Critical density factors
	const CritOverDensity = 577                 # For cosmological simulations
	const CritPhysDensity = 0.318               # For newtonian simulations
	
end;

# ╔═╡ b4f80502-cb91-4358-9229-414ad8774ded
begin
	
	# T [internal_units] * t_s = T [s]
	t_s = UnitTime_in_s / HubbleParam
	
	# T [internal_units] * t_Gyr = T [Gyr]
	t_Gyr = t_s / (SEC_PER_YEAR * 1e9)
	
	# RHO [internal_units] * rho_cgs = RHO [g cm^(-3)]
    rho_cgs = UnitDensity_in_cgs * HubbleParam * HubbleParam * a3inv
	
	# RHO [internal_units] * rho_cosmo = RHO [Mₒ pc^(-3)]
    rho_cosmo = rho_cgs * PARSEC * PARSEC * PARSEC / SOLAR_MASS
	
end;

# ╔═╡ f492e1c3-b2e8-4dc3-8a7c-fd4499f4ef04
1.0u"G" |> u"cm^3/g/s^2"

# ╔═╡ 4bbac9f0-81ad-11eb-3bf6-d1280cead58d
md"### Marios's $R_H$ from $\tau_C = 797\,\text{Myr}$"

# ╔═╡ aecf0420-81ad-11eb-2fd2-438f479be225
((u"k" * 100u"K") / (797UnitfulAstro.Myr * π * u"G" * UnitfulAstro.Msun^2 * UnitfulAstro.pc^(-4))) |> u"cm^3/s"

# ╔═╡ 33a5723e-81b0-11eb-17c3-4b445936501d
md"### Evaluation of my $C_1$, which is equal to Carlos' $\tau_w$"

# ╔═╡ fe312110-81ad-11eb-0be5-d92e515464be
(4 * u"k" * Tw) / ((4.1e-10 * 1e4^(-0.8))u"cm^3/s" * π * u"G" * UnitfulAstro.Msun^2 * UnitfulAstro.pc^(-4)) |> UnitfulAstro.Myr

# ╔═╡ 4ce8e24e-81b0-11eb-18ea-01270d41ec2c
md"### Recombination rate, Cecilia $[\text{cm}^3 \, \text{s}^{-1}]$"

# ╔═╡ 7e5c4180-81ae-11eb-335b-1fba28a5f8d7
f(1e4)

# ╔═╡ 5be99560-81c4-11eb-18c9-a37897343dcc
md"### Recombination rate, for the specific value $T = T_\text{warm}$" 

# ╔═╡ 815165a0-81ae-11eb-3a29-8504d93561d9
f(1e4) |> u"cm^3" / UnitfulAstro.Gyr

# ╔═╡ 0bcbeb40-81b1-11eb-0e99-cd107a40dca6
md"### $C_2$ constant"

# ╔═╡ bc161af0-81ae-11eb-1827-6f9a6e55487d
((4 * u"k" * Tw) / (π * u"G" * f(1e4))) |> UnitfulAstro.Msun^2 * UnitfulAstro.pc^(-4) * UnitfulAstro.Gyr

# ╔═╡ 41590680-81b1-11eb-33d3-b10a16e409a4
md"### $C_2$ constant $\rightarrow$ Carlos' $\tau_w$"

# ╔═╡ d0d8e4e0-81ae-11eb-2de4-578fafcd3b80
((4 * u"k" * Tw) / (π * u"G" * f(1e4))) |> UnitfulAstro.Msun^2 * UnitfulAstro.pc^(-4) * UnitfulAstro.Myr

# ╔═╡ 62e26210-81b1-11eb-1194-cdaaf93b135e
md"### $R_h$ constant in better units"

# ╔═╡ dedf97a0-81ae-11eb-193b-d3d2023c2fb0
Rh = 6e-17u"cm^3/s" |> u"cm^3" / UnitfulAstro.Gyr

# ╔═╡ b9aeb082-81b1-11eb-2c62-1ba2b4bb2087
md"### $C_4$ constant with $Z_\odot$ = 0.0134"

# ╔═╡ ecf7d690-81ae-11eb-0f3d-63b9191d8843
((0.0134 * u"k" * Tc) / (π * u"G" * Rh)) |> UnitfulAstro.Msun^2 * UnitfulAstro.pc^(-4) * UnitfulAstro.Gyr

# ╔═╡ fc8d8660-81b1-11eb-0bbd-13ce824a4ac5
md"### $C_4$ constant with $Z_\odot$ = 0.006 (Carlos)"

# ╔═╡ f306bdd2-81ae-11eb-3b31-d7052a1397ab
((0.006 * u"k" * Tc) / (π * u"G" * Rh)) |> UnitfulAstro.Msun^2 * UnitfulAstro.pc^(-4) * UnitfulAstro.Gyr

# ╔═╡ 563a115e-81b2-11eb-2832-9325d67f3880
md"### $C_4$ constant without $Z_\odot$"

# ╔═╡ fabe4c00-81ae-11eb-09e0-517f857d9781
((u"k" * Tc) / (π * u"G" * Rh)) |> UnitfulAstro.Msun^2 * UnitfulAstro.pc^(-4) * UnitfulAstro.Gyr

# ╔═╡ 6fc8460e-81b2-11eb-104b-47257a3a604e
md"### $C_4$ constant $\rightarrow$ Carlos' $\tau_C$"

# ╔═╡ 02cc1080-81af-11eb-1cff-03bad4d1dbba
((u"k" * Tc) / (π * u"G" * Rh)) |> UnitfulAstro.Msun^2 * UnitfulAstro.pc^(-4) * UnitfulAstro.Myr

# ╔═╡ 8b9a2930-81b2-11eb-3e46-8d619500b9f8
md"### massFactor - From proton density to mass area density"

# ╔═╡ 1e67ef80-81af-11eb-2e44-e3f1ebe11780
massFactor = galaxythickness * u"mp/cm^3" |> UnitfulAstro.Msun / UnitfulAstro.pc^2

# ╔═╡ 98585660-81b2-11eb-2df9-2bbcab54c352
md"### massFactor2 - From proton density to mass area density"

# ╔═╡ 3b049cb0-81af-11eb-3438-8154b464c2e1
u"mp/cm^3" |> UnitfulAstro.Msun / UnitfulAstro.pc^3

# ╔═╡ d53abb90-81b2-11eb-221f-872821b1ac19
md"### $C_1$ and $C_2$ for the Ascasibar et al. model"

# ╔═╡ ecd42522-81b2-11eb-1543-8dca5a15b266
md"Constant $C_1$"

# ╔═╡ 0f142de0-81b0-11eb-2c95-75ed772278bc
((2 * u"k" * (Tw + Teff)) / (π * u"G" * σνB)) |> UnitfulAstro.Msun^2 * UnitfulAstro.pc^(-4) * UnitfulAstro.Gyr

# ╔═╡ 0a4930f2-81b3-11eb-2913-c19ffcb4ae00
md"Constant $C_2$"

# ╔═╡ 1c480040-81b0-11eb-2e4f-87b69a6b813e
((Zsun * u"k" * *(Tc + Teff)) / (3 * π * Rh * u"G")) |> UnitfulAstro.Msun^2 * UnitfulAstro.pc^(-4) * UnitfulAstro.Gyr

# ╔═╡ ace59bf0-81b3-11eb-3875-614fb6a2f56d
md"### $K_1$ constant"

# ╔═╡ bd619470-81b3-11eb-3733-05f5c0229b52
K1 = (π * u"G") / (2 * u"k" * (Tw + Teff)) |> UnitfulAstro.pc^4 * UnitfulAstro.Msun^(-2) * u"cm^-3" 

# ╔═╡ 00af58c0-81b4-11eb-24a3-b9d0b8a00829
md"### $K_2$ constant"

# ╔═╡ 04a30120-81b4-11eb-32a3-6346f5a572a1
K2 = (3 * π * u"G") / (2 * u"k" * (Tc + Teff)) |> UnitfulAstro.pc^4 * UnitfulAstro.Msun^(-2) * u"cm^-3" 

# ╔═╡ 1a2e4b30-81b4-11eb-1973-3dd9d03f8ec8
md"### $g_0$ for the comparison between the basic model and the Ascasibar et al. model (high density)"

# ╔═╡ 3bf265d0-81b4-11eb-1acb-310d0bba54dd
begin
	local n = 1000u"cm^-3"
	sqrt(2 * n / (K1 + K2)) |> UnitfulAstro.Msun * UnitfulAstro.pc^-2
end

# ╔═╡ 9869a990-81b4-11eb-3ee9-1b839a36fccf
md"### $g_0$ for the comparison between the basic model and the Ascasibar et al. model (low density)"

# ╔═╡ a723589e-81b4-11eb-1a8d-1553977b5b2e
begin
	local n = 1u"cm^-3"
	sqrt(2 * n / (K1 + K2)) |> UnitfulAstro.Msun * UnitfulAstro.pc^-2
end

# ╔═╡ 6cfd10c0-81b5-11eb-120d-058a16e80989
md"### $P_1$ constant (conversion from mass density to number density, $n_e = n_i = P_1\,i(t)$, where $[i(t)] = \text{M}_\odot \, \text{pc}^{-3}$)"

# ╔═╡ aefbe6e0-81b5-11eb-06df-65cd61a32c66
UnitfulAstro.Msun * UnitfulAstro.pc^-3 * u"mp^-1" |> u"cm^-3"

# ╔═╡ d6fb28e0-81b5-11eb-02ff-4f280e4e0cb2
md"### $P_2$ and $P_3$ constant (conversion from mass density to number density, $n_H = n_a+ 2 \, n_m = P_2 \, a(t) + P_3 \, 2 \, m(t)$, where $[a(t)] = [m(t)] = \text{M}_\odot \, \text{pc}^{-3}$)"

# ╔═╡ 509033e0-81ba-11eb-1da9-590737aea01f
mH = u"molH" * u"mol^-1*Na^-1" |> u"g"

# ╔═╡ 2091a922-81b6-11eb-1fc9-5ba004c2751a
UnitfulAstro.Msun * UnitfulAstro.pc^-3 * mH^-1 |> u"cm^-3"

# ╔═╡ baef6252-81bb-11eb-06ad-0d930a782a02
 @compound H2

# ╔═╡ b72ba38e-81bb-11eb-0d0c-e167d532dc79
mH2 = 1molH₂ * u"mol^-1*Na^-1" |> u"g"

# ╔═╡ f408ffb0-81bb-11eb-2058-fbdc0dd253a0
UnitfulAstro.Msun * UnitfulAstro.pc^-3 * mH2^-1 |> u"cm^-3"

# ╔═╡ 844fb7e0-81b6-11eb-39f2-61e8c847db1e
md"### $C_1$ and $C_2$ for the Self consistent Ascasibar et al. model"

# ╔═╡ 0f587e50-81ba-11eb-00f4-b140f5717e97
md"Constant $C_1$"

# ╔═╡ 1ad44070-81ba-11eb-0662-bb47b0103d7d
P1 = UnitfulAstro.Msun * UnitfulAstro.pc^-3 * u"mp^-1" |> u"cm^-3"

# ╔═╡ 14440950-81bc-11eb-13a3-1790199caaf6
1 / (P1 * σνB) |> UnitfulAstro.Gyr    # C1

# ╔═╡ 31d88400-81bc-11eb-024f-9704e887fbfe
md"Constant $C_2$"

# ╔═╡ fff759c0-81bb-11eb-0225-0140ad03b022
P2 = UnitfulAstro.Msun * UnitfulAstro.pc^-3 * mH^-1 |> u"cm^-3"

# ╔═╡ 3b8b1ee0-81bc-11eb-0b13-a34201b0187b
Zsun / (2 * Rh * P2) |> UnitfulAstro.Gyr    # C2

# ╔═╡ 5b249150-81bc-11eb-165d-f58f952906f8
md"### HUBBLE_CONST for GADGETPlotting.jl"

# ╔═╡ e334b890-81bc-11eb-2fd0-0b6e15e93fc3
100u"km*s^-1" * UnitfulAstro.Mpc^-1 |> UnitfulAstro.Gyr^-1

# ╔═╡ 2d5b1cc0-81bd-11eb-22b3-03d73dd827e4
md"### Constant in $\tau_S$ for the fully coupled model with variable $\tau_S$"

# ╔═╡ 41f09520-81bd-11eb-2e67-4de804b697f1
begin
	local eff = 0.01
	(1 / eff) * sqrt(3 * π / 32u"G") |> UnitfulAstro.Gyr * UnitfulAstro.Msun^(1/2) * UnitfulAstro.pc^(-3/2)
end

# ╔═╡ 4685db10-81c5-11eb-3324-7d72858f3647
md"### Unit conversion constants"

# ╔═╡ 54e8b792-81c5-11eb-2f8c-db1704e37605
tCGS = 9.8e8UnitfulAstro.yr |> u"s"

# ╔═╡ 724a3b62-81c5-11eb-0db6-79b6a25af7bb
tGyr = 9.8e8UnitfulAstro.yr |> UnitfulAstro.Gyr

# ╔═╡ 927d7cd2-81c5-11eb-1661-6f106df93c07
ρCGS = 1e10 * UnitfulAstro.Msun * UnitfulAstro.kpc^-3 |> u"g*cm^-3"

# ╔═╡ d1fcfb10-81c5-11eb-1925-6ddbca209dfb
ρCosmo = 1e10 * UnitfulAstro.Msun * UnitfulAstro.kpc^-3 |> UnitfulAstro.Msun * UnitfulAstro.pc^-3

# ╔═╡ b97b4980-8598-11eb-1aaa-2dff6b61d37c
md"### Conversions from internal units"

# ╔═╡ cfe5a8ee-8598-11eb-3958-81452b5ce8db
10 * rho_cosmo # [Mₒ pc^(-3)]

# ╔═╡ 5bd8ae69-a89d-472d-899a-1b73e5cd82a9
10 * rho_cosmo * P1 # [cm^(-3)]

# ╔═╡ 8300f32a-fb20-46ea-9cec-234cb30fdc55
10 * rho_cosmo * P1 * massFactor # [Mₒ pc^(-2)]

# ╔═╡ 7195d6e2-1ac2-426a-b307-54748754faa7
md"### Particle number to mass density"

# ╔═╡ 669b36f4-c18b-4624-8e9a-71bd15b0f0d2
1u"mp/cm^3" |> UnitfulAstro.Msun / UnitfulAstro.kpc^3

# ╔═╡ 2cc04d53-d45a-42e9-af10-5ad186db6a06
10u"mp/cm^3" |> UnitfulAstro.Msun / UnitfulAstro.kpc^3

# ╔═╡ 2b92d2d8-e1ee-4ca6-9d5b-7dffdc685120
100u"mp/cm^3" |> UnitfulAstro.Msun / UnitfulAstro.kpc^3

# ╔═╡ ab8bf3fa-581a-483b-9945-f3e94725bc6e
1000u"mp/cm^3" |> UnitfulAstro.Msun / UnitfulAstro.kpc^3

# ╔═╡ 2b113886-797e-4e21-b2ef-867ebb7e21ab
5000u"mp/cm^3" |> UnitfulAstro.Msun / UnitfulAstro.kpc^3

# ╔═╡ d9d2f2ee-588e-4f0f-9e9f-6db835121ab8
10000u"mp/cm^3" |> UnitfulAstro.Msun / UnitfulAstro.kpc^3

# ╔═╡ 428a6794-8f79-4166-b1c9-4eaa8dda1f61
md"### Mean particle distance and smoothing"

# ╔═╡ 915d0b51-3a8d-4ccd-94e7-e6165edd2254
function mean_distance(R::Float64, packing_eff::Float64, l_n::Int64)
	N = l_n * l_n * l_n
	r = R * (packing_eff / N)^(1 / 3.0)
	
	return 2 * r
	
end;

# ╔═╡ 516186c2-8bb6-4ea6-9f28-2e25622110c0
md"__Irregular packing__"

# ╔═╡ 06b6fb35-69c9-4ffb-97c8-a3159e0b57c8
mean_distance(100.0, 0.64, 32) / 20

# ╔═╡ abc4e9e0-1da4-44ff-b599-9a59bfbf5d5d
mean_distance(100.0, 0.64, 32) / 50

# ╔═╡ bee89c62-06f1-4711-9c6c-5e3a696cf301
md"__Face-centred cubic (FCC)__"

# ╔═╡ 40886e8c-527f-4057-8635-15ff2780f192
mean_distance(100.0, 0.74048, 32) / 20

# ╔═╡ 9a24c5e3-d92a-4fa5-9a13-3bebd60900e0
mean_distance(100.0, 0.74048, 32) / 50

# ╔═╡ 176f6707-4327-446d-b36f-0f328ad97bec
md"__Cubic lattice__"

# ╔═╡ aec3a280-228b-4f6d-aced-e2c26bcd1eb1
mean_distance(100.0, 0.5236, 32) / 20

# ╔═╡ b78f3224-d6dc-4ea2-98aa-b68e993be652
mean_distance(100.0, 0.5236, 32) / 50

# ╔═╡ 8cd37573-1b63-4fbb-b813-3caf07368a36
md"__For diferent resolutions__"

# ╔═╡ ea217cc6-61a3-46e8-88eb-361dd2cb07bf
factor_64 = mean_distance(100.0, 0.64, 64) / mean_distance(100.0, 0.64, 32)

# ╔═╡ 3b498ba3-f5e3-4e61-b9d3-761db7cc3ee3
factor_128 = mean_distance(100.0, 0.64, 128) / mean_distance(100.0, 0.64, 32)

# ╔═╡ 75faa6bb-77cc-4372-bb98-3c3c82dfc481
md"### Critical density"

# ╔═╡ 7795182f-a99e-443b-b579-cf0e4ff464bf
PhysDensThresh = (CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / UnitDensity_in_cgs) * 1.0e10UnitfulAstro.Msun / UnitfulAstro.kpc^3

# ╔═╡ 65f5eb7c-484c-4cc1-ae64-09d321180a9e
(CritPhysDensity / HYDROGEN_MASSFRAC)u"cm^-3"

# ╔═╡ 05931e68-2f9d-41d8-9032-e21f7eff0082
md"CritPhysDensity = 0.150"

# ╔═╡ 752022ab-04f4-49a3-bc0f-fed555889241
(0.150 * PROTONMASS / HYDROGEN_MASSFRAC / UnitDensity_in_cgs) * 1.0e10UnitfulAstro.Msun / UnitfulAstro.kpc^3

# ╔═╡ 17f8cb5a-7c94-4d79-8f48-1630e92ec02c
(0.150 / HYDROGEN_MASSFRAC)u"cm^-3"

# ╔═╡ e99a22a6-adc5-49c6-a826-77eb1b4fa6dd
md"CritPhysDensity = 0.0.05"

# ╔═╡ e7ca6285-447b-4686-83bb-5cf3ba4aed41
(0.05 * PROTONMASS / HYDROGEN_MASSFRAC / UnitDensity_in_cgs) * 1.0e10UnitfulAstro.Msun / UnitfulAstro.kpc^3

# ╔═╡ 83091cb6-dc66-4027-8af4-54bffd95ae9d
(0.05 / HYDROGEN_MASSFRAC)u"cm^-3"

# ╔═╡ bed5dcad-8f55-4a3f-aafc-471690338eb4
md"Cosmological case"

# ╔═╡ 7d9c0a99-c5df-499a-b852-7d13802c2b73
OverDensThresh = CritOverDensity * OmegaBaryon * 3 * HUBBLE^2 / (8π * GRAVITY)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
UnitfulMoles = "999f2bd7-36bf-5ba7-9bc1-c9473aa75374"

[compat]
PlutoUI = "~0.7.9"
Unitful = "~1.9.0"
UnitfulAstro = "~1.0.1"
UnitfulMoles = "~0.1.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a981a8ef8714cba2fd9780b22fd7a469e7aaf56d"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.9.0"

[[UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "dd21b5420bf6e9b76a8c6e56fb575319e7b1f895"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.6.1"

[[UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "b154be4ca9610e9c9dbb9dba98b2bd750539630b"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.0.1"

[[UnitfulMoles]]
deps = ["Unitful"]
git-tree-sha1 = "0ce82f066f57a9c2dd00ee5c47c9a3677be61e70"
uuid = "999f2bd7-36bf-5ba7-9bc1-c9473aa75374"
version = "0.1.0"
"""

# ╔═╡ Cell order:
# ╟─b8930ef2-e0a3-4ff9-9cce-88dcaf486cef
# ╠═8ae23d20-81ad-11eb-0644-4796c8c1e8c0
# ╟─bf6ad952-81b0-11eb-1de7-57344a180233
# ╠═44bcd24e-81ae-11eb-392f-43c49f59f6da
# ╠═21d3bdb0-81c4-11eb-3e58-973170b3ce20
# ╟─f68db520-81c3-11eb-0463-412bda3df6c6
# ╠═47281080-81af-11eb-0bc3-8762ecd9a5f0
# ╠═b4f80502-cb91-4358-9229-414ad8774ded
# ╠═f492e1c3-b2e8-4dc3-8a7c-fd4499f4ef04
# ╟─4bbac9f0-81ad-11eb-3bf6-d1280cead58d
# ╠═aecf0420-81ad-11eb-2fd2-438f479be225
# ╟─33a5723e-81b0-11eb-17c3-4b445936501d
# ╠═fe312110-81ad-11eb-0be5-d92e515464be
# ╟─4ce8e24e-81b0-11eb-18ea-01270d41ec2c
# ╠═7e5c4180-81ae-11eb-335b-1fba28a5f8d7
# ╟─5be99560-81c4-11eb-18c9-a37897343dcc
# ╠═815165a0-81ae-11eb-3a29-8504d93561d9
# ╟─0bcbeb40-81b1-11eb-0e99-cd107a40dca6
# ╠═bc161af0-81ae-11eb-1827-6f9a6e55487d
# ╟─41590680-81b1-11eb-33d3-b10a16e409a4
# ╠═d0d8e4e0-81ae-11eb-2de4-578fafcd3b80
# ╟─62e26210-81b1-11eb-1194-cdaaf93b135e
# ╠═dedf97a0-81ae-11eb-193b-d3d2023c2fb0
# ╟─b9aeb082-81b1-11eb-2c62-1ba2b4bb2087
# ╠═ecf7d690-81ae-11eb-0f3d-63b9191d8843
# ╟─fc8d8660-81b1-11eb-0bbd-13ce824a4ac5
# ╠═f306bdd2-81ae-11eb-3b31-d7052a1397ab
# ╟─563a115e-81b2-11eb-2832-9325d67f3880
# ╠═fabe4c00-81ae-11eb-09e0-517f857d9781
# ╟─6fc8460e-81b2-11eb-104b-47257a3a604e
# ╠═02cc1080-81af-11eb-1cff-03bad4d1dbba
# ╟─8b9a2930-81b2-11eb-3e46-8d619500b9f8
# ╠═1e67ef80-81af-11eb-2e44-e3f1ebe11780
# ╟─98585660-81b2-11eb-2df9-2bbcab54c352
# ╠═3b049cb0-81af-11eb-3438-8154b464c2e1
# ╟─d53abb90-81b2-11eb-221f-872821b1ac19
# ╟─ecd42522-81b2-11eb-1543-8dca5a15b266
# ╠═0f142de0-81b0-11eb-2c95-75ed772278bc
# ╟─0a4930f2-81b3-11eb-2913-c19ffcb4ae00
# ╠═1c480040-81b0-11eb-2e4f-87b69a6b813e
# ╟─ace59bf0-81b3-11eb-3875-614fb6a2f56d
# ╠═bd619470-81b3-11eb-3733-05f5c0229b52
# ╟─00af58c0-81b4-11eb-24a3-b9d0b8a00829
# ╠═04a30120-81b4-11eb-32a3-6346f5a572a1
# ╟─1a2e4b30-81b4-11eb-1973-3dd9d03f8ec8
# ╠═3bf265d0-81b4-11eb-1acb-310d0bba54dd
# ╟─9869a990-81b4-11eb-3ee9-1b839a36fccf
# ╠═a723589e-81b4-11eb-1a8d-1553977b5b2e
# ╟─6cfd10c0-81b5-11eb-120d-058a16e80989
# ╠═aefbe6e0-81b5-11eb-06df-65cd61a32c66
# ╟─d6fb28e0-81b5-11eb-02ff-4f280e4e0cb2
# ╟─509033e0-81ba-11eb-1da9-590737aea01f
# ╠═2091a922-81b6-11eb-1fc9-5ba004c2751a
# ╟─baef6252-81bb-11eb-06ad-0d930a782a02
# ╟─b72ba38e-81bb-11eb-0d0c-e167d532dc79
# ╠═f408ffb0-81bb-11eb-2058-fbdc0dd253a0
# ╟─844fb7e0-81b6-11eb-39f2-61e8c847db1e
# ╟─0f587e50-81ba-11eb-00f4-b140f5717e97
# ╠═1ad44070-81ba-11eb-0662-bb47b0103d7d
# ╠═14440950-81bc-11eb-13a3-1790199caaf6
# ╟─31d88400-81bc-11eb-024f-9704e887fbfe
# ╠═fff759c0-81bb-11eb-0225-0140ad03b022
# ╠═3b8b1ee0-81bc-11eb-0b13-a34201b0187b
# ╟─5b249150-81bc-11eb-165d-f58f952906f8
# ╠═e334b890-81bc-11eb-2fd0-0b6e15e93fc3
# ╟─2d5b1cc0-81bd-11eb-22b3-03d73dd827e4
# ╠═41f09520-81bd-11eb-2e67-4de804b697f1
# ╟─4685db10-81c5-11eb-3324-7d72858f3647
# ╠═54e8b792-81c5-11eb-2f8c-db1704e37605
# ╠═724a3b62-81c5-11eb-0db6-79b6a25af7bb
# ╠═927d7cd2-81c5-11eb-1661-6f106df93c07
# ╠═d1fcfb10-81c5-11eb-1925-6ddbca209dfb
# ╠═b97b4980-8598-11eb-1aaa-2dff6b61d37c
# ╠═cfe5a8ee-8598-11eb-3958-81452b5ce8db
# ╠═5bd8ae69-a89d-472d-899a-1b73e5cd82a9
# ╠═8300f32a-fb20-46ea-9cec-234cb30fdc55
# ╟─7195d6e2-1ac2-426a-b307-54748754faa7
# ╠═669b36f4-c18b-4624-8e9a-71bd15b0f0d2
# ╠═2cc04d53-d45a-42e9-af10-5ad186db6a06
# ╠═2b92d2d8-e1ee-4ca6-9d5b-7dffdc685120
# ╠═ab8bf3fa-581a-483b-9945-f3e94725bc6e
# ╠═2b113886-797e-4e21-b2ef-867ebb7e21ab
# ╠═d9d2f2ee-588e-4f0f-9e9f-6db835121ab8
# ╟─428a6794-8f79-4166-b1c9-4eaa8dda1f61
# ╠═915d0b51-3a8d-4ccd-94e7-e6165edd2254
# ╟─516186c2-8bb6-4ea6-9f28-2e25622110c0
# ╠═06b6fb35-69c9-4ffb-97c8-a3159e0b57c8
# ╠═abc4e9e0-1da4-44ff-b599-9a59bfbf5d5d
# ╟─bee89c62-06f1-4711-9c6c-5e3a696cf301
# ╠═40886e8c-527f-4057-8635-15ff2780f192
# ╠═9a24c5e3-d92a-4fa5-9a13-3bebd60900e0
# ╟─176f6707-4327-446d-b36f-0f328ad97bec
# ╠═aec3a280-228b-4f6d-aced-e2c26bcd1eb1
# ╠═b78f3224-d6dc-4ea2-98aa-b68e993be652
# ╟─8cd37573-1b63-4fbb-b813-3caf07368a36
# ╠═ea217cc6-61a3-46e8-88eb-361dd2cb07bf
# ╠═3b498ba3-f5e3-4e61-b9d3-761db7cc3ee3
# ╟─75faa6bb-77cc-4372-bb98-3c3c82dfc481
# ╠═7795182f-a99e-443b-b579-cf0e4ff464bf
# ╠═65f5eb7c-484c-4cc1-ae64-09d321180a9e
# ╠═05931e68-2f9d-41d8-9032-e21f7eff0082
# ╠═752022ab-04f4-49a3-bc0f-fed555889241
# ╠═17f8cb5a-7c94-4d79-8f48-1630e92ec02c
# ╟─e99a22a6-adc5-49c6-a826-77eb1b4fa6dd
# ╠═e7ca6285-447b-4686-83bb-5cf3ba4aed41
# ╠═83091cb6-dc66-4027-8af4-54bffd95ae9d
# ╟─bed5dcad-8f55-4a3f-aafc-471690338eb4
# ╠═7d9c0a99-c5df-499a-b852-7d13802c2b73
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

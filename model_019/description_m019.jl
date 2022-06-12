### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ 982068e0-59bb-11ec-27f5-51126c2ba1df
using DifferentialEquations, PlutoUI, TikzPictures, Unitful, UnitfulAstro, SpecialFunctions, QuadGK, CairoMakie, LaTeXStrings, DataFrames, DataFramesMeta, DelimitedFiles, Trapz, LsqFit, Interpolations, CSV, HDF5, Symbolics, Statistics

# ╔═╡ 08df960b-fd82-43ba-a9dc-bf5e83af587e
TableOfContents(title="🌌 Model 019", depth=4)

# ╔═╡ 3b1726c6-60c2-45be-932d-efa8d2ef23e0
Markdown.MD(
	Markdown.Admonition(
		"note", 
		"Names", 
		[md"Model 019"]
	)
)

# ╔═╡ 6173713c-92b0-43ec-8713-1cbf442aa1ce
md"# Description"

# ╔═╡ a842b24e-8d26-41ab-9de3-91632aede893
md"""

## Phases and notation

Because the interstellar medium (ISM) in a galaxy is not homogeneous, we will model it as a multiphase structure. For that purpose, we consider an ionized phase with a temperature of $\sim \! 10^4 \, \mathrm{K}$; an atomic component whose temperature is $\sim \! 100 \, \mathrm{K}$; and a molecular phase with a temperature of $\sim \! 10 \, \mathrm{K}$. The other two phases considered are the stars and all other elements present in the system (metals).
Metals interact only with the stars and do not contribute to mass balance. 

For every reaction, we will take into account only the dominant channel as a first approximation. So, even though the gas contains Hydrogen and Helium, only Hydrogen reactions are incorporated into the model.

The phases interchange mass with each other during the evolution of the SPH particle under study. The processes involved are the photoionization of atoms; the recombination of electrons with ions; the conversion of atomic hydrogen into molecular hydrogen; and the destruction of the latter owing to the photodissociation caused by ultraviolet (UV) light from the young stellar population. In addition, we consider the formation by supernovas of metals and ionized gas, and the influence of molecular and atomic gas on the star formation rate.

We characterized the mass of the phases by their density fraction, namely

*  Ionized gas $\rightarrow \ i_f(t) ≔ \rho_i / \rho \, ,$
*  Atomic gas $\rightarrow \ a_f(t) ≔ \rho_a / \rho \, ,$
*  Molecular gas $\rightarrow \ m_f(t) ≔ \rho_m / \rho \, ,$
*  Stars $\rightarrow \ s_f(t) ≔ \rho_s / \rho \, ,$
*  Metals $\rightarrow \ z_f(t) ≔ \rho_z / \rho \, ,$

where $\rho_i$, $\rho_a$, $\rho_m$, $\rho_s$ and $\rho_z$ are the volume densities of each phase and

$\begin{align}
	\rho_g(t) &= \rho_i(t) + \rho_a(t) + \rho_m(t) \, , \\
	g_f(t) &= \rho_g(t) / \rho \, , \\
    \rho(t) &= \rho_g(t) + \rho_s(t) \, ,
\end{align}$

are the gas density, the gas density fraction, and the total density, of the SPH particle in consideration.

Notice that $z_f$ is not the metallicity as is commonly defined in the literature, where in general $z = \rho_z / \rho_g$. But, given that the stellar fraction (in the short times that are relevant for us) will always be very small, the difference is negligible.
"""

# ╔═╡ 64787011-b5b8-42be-b6e4-37ebc5138b3e
md"## Relations between the phases"

# ╔═╡ 76fe97bd-36c8-40d2-9b5a-0ea5059bd7c7
md"""
The following diagram shows the relation between the different phases,
"""

# ╔═╡ 14c7f574-0623-4254-b8f7-97984d32351c
TikzPicture(
	L"""
		\node[box, white] (stars) {Stars};
		\node[box, white, text width=2em, above=2cm of stars] (atom) {H};
		\node[box, white, text width=2em, right=2cm of atom] (molecule) {\ch{H2}};
		\node[box, white, text width=2em, left=2cm of atom] (ion) {\ch{H+}};
		\node[box, white, below=1cm of stars] (metals) {Metals};
		\draw[line, white, ->]
		(ion) edge [bend left, "\textcolor{red}{recombination}"] (atom)
		(atom) edge [bend left, "\textcolor{blue}{condensation}"] (molecule)
		(molecule) edge [bend left,"\textcolor{new_green}{dissociation}"] (atom)
		(atom) edge [bend left,"\textcolor{violet}{ionization}"] (ion)
		(stars) edge [bend left, "\textcolor{cyan}{supernova}"] (ion)
		(molecule) edge [bend left, "\textcolor{yellow}{star formation}"] (stars)
		(stars.west) edge [bend right] node [midway, left] {\textcolor{pink}{\small{supernova}}} (metals.west)
		(metals.east) to [bend right] node [midway, right] {\textcolor{orange}{\small{star formation}}} (stars.east);
	""", 
	width="55em",
	preamble = """
		\\usepackage{chemformula}
		\\definecolor{new_green}{HTML}{00AC06}
		\\definecolor{yellow}{HTML}{B59E15}
		\\definecolor{violet}{HTML}{7F13EC}
		\\definecolor{cyan}{HTML}{0BC3DF}
		\\definecolor{pink}{HTML}{FF0B8C}
		\\definecolor{orange}{HTML}{F39C12}
		\\usetikzlibrary{shapes.misc, arrows, positioning, quotes, fit}
		\\tikzset{
    		>=stealth',
    		box/.style={
        		rectangle,
        		rounded corners,
        		draw=black, 
        		very thick,
        		text width=4em,
        		minimum height=2em,
        		text centered,
    		},
			line/.style = {
				thick,
			},
			every edge quotes/.append style = {
				font=\\small, 
				align=center, 
				auto,
			},
			myrect/.style={
				rectangle,
				draw,
				inner sep=0pt,
				fit=#1,
				thick, 
				rounded corners,
			},
		}
	""",
)

# ╔═╡ bc9ab101-7cc3-4baa-b83d-ce546f6b576d
md"## Equations"

# ╔═╡ 047bbd39-9cf9-4bd7-b38e-16aa505b0b08
md"""
Given the previous interaction diagram we have a system of five first-order ODEs:
"""

# ╔═╡ 2fe0dc4c-da44-4fc8-bef8-1fa615a0fe4a
TikzPicture(
	L"""
	\node[white] {
  	${\boldmath\begin{aligned}
		\dv{}{t}i_f(t) &= - \textcolor{red}{\frac{i_f(t)}{\tau_R(t)}} + \textcolor{violet}{\eta_\text{ion}\,\psi(t)} + \textcolor{cyan}{R\,\psi(t)} \, , \\
	\dv{}{t}a_f(t) &= - \textcolor{blue}{\frac{a_f(t)}{\tau_C(t)}} + \textcolor{red}{\frac{i_f(t)}{\tau_R(t)}} + \textcolor{new_green}{\eta_\text{diss}\,\psi(t)} - \textcolor{violet}{\eta_\text{ion}\,\psi(t)} \, , \\
	\dv{}{t}m_f(t) &= \textcolor{blue}{\frac{a_f(t)}{\tau_C(t)}} - \textcolor{new_green}{\eta_\text{diss}\,\psi(t)} - \textcolor{yellow}{\psi(t)} \, , \\
	\dv{}{t}s_f(t) &= \textcolor{yellow}{\psi(t)} - \textcolor{cyan}{R\,\psi(t)} \, , \\
	\dv{}{t}z_f(t) &= \textcolor{pink}{Z_{SN}\,R\,\psi(t)} - \textcolor{orange}{z_f\,\psi(t)} \, ,
	\end{aligned}}$};
	""", 
	width="35em",
	preamble = """
		\\usepackage{chemformula}
		\\usepackage{physics}
		\\setlength{\\jot}{10pt}
		\\definecolor{new_green}{HTML}{00AC06}
		\\definecolor{yellow}{HTML}{B59E15}
		\\definecolor{violet}{HTML}{7F13EC}
		\\definecolor{cyan}{HTML}{0BC3DF}
		\\definecolor{pink}{HTML}{FF0B8C}
		\\definecolor{orange}{HTML}{F39C12}
	""",
)

# ╔═╡ eaf272c7-4162-4a9a-92e3-9835c6158394
md"""
where

*  $\tau_R$ is the time scale of atomic gas formation from ionized gas, generally called recombination time.

*  $\tau_C$ is the time scale of molecular gas formation from atomic gas, generally called condensation (or cloud formation) time.

*  $\psi(t)$ is the star formation rate (SFR).

*  $\eta_\text{ion}$ is the rate of atomic gas ionization by stars, per unit of created stellar mass.

*  $\eta_\text{diss}$ is the rate of molecular gas dissociation by stars, per unit of created stellar mass.

*  $R$ is the mass of ionized gas produced per unit of created stellar mass. The model assumes instantaneous death for stars with more than eight solar masses.

*  $Z_{SN}$ is the mass of metals produced per unit of ionized gas created during the stellar life cycle. Then, $Z_{SN}\,R$ would be the mass of metals produced per unit of created stellar mass.

From the previous equations we explicitly see mass conservation for $i_f$, $a_f$, $m_f$ and $s_f$,

$\begin{equation}
    \sum_j \frac{\mathrm{d}}{\mathrm{d}t}j(t) = 0 \, , \qquad j = i_f, a_f, m_f, s_f \, .
\end{equation}$
"""

# ╔═╡ ac553b12-4857-4cc1-8ea2-fe9e8863b429
md"""

## Star formation rate

For the star formation rate we take into account the strong correlation between molecular Hydrogen and star formation ([Bigiel2008](https://doi.org/10.1088/0004-6256/136/6/2846), [Bigiel2010](https://doi.org/10.1088/0004-6256/140/5/1194), [Wong2002](https://doi.org/10.1086/339287), [Robertson2008](https://doi.org/10.1086/587796), [Halle2013](https://doi.org/10.1051/0004-6361/201220952), [Thompson2013](https://doi.org/10.1088/0004-637x/780/2/145)). In particular we will follow [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) with

$\begin{equation}
	\psi = \frac{\alpha \, a_f + \beta \, m_f}{\tau_S}
\end{equation}$
where $\alpha$ y $\beta$ are dimensionless free parameters.

With the choice $\alpha = 0$ y $\beta = 1$ we get

$\begin{equation}
	\psi = \frac{m_f}{\tau_S}
\end{equation}$

where $\tau_S$ is the characteristic timescale for star formation.
"""

# ╔═╡ dc6fd12b-c821-4e20-a896-25c8aab9df94
md"## Time parameters"

# ╔═╡ 1d27ec35-65ca-4c94-9e8d-54d1c11e759f
md"""

### Star formation time

Following [Krumholz2019](https://doi.org/10.1146/annurev-astro-091918-104430), we define $\tau_S$ as the characteristic timescale for star formation; i.e., all the gas which can be converted to stars has done so in a timescale $\tau_S$. In particular, we have

$\begin{equation}
    \tau_S = \frac{\epsilon_\star}{\epsilon_\text{ff}}\,t_\text{ff} \, ,
\end{equation}$
where $\epsilon_\star$ is the fraction of total gas mass which will be converted to stars, $t_\text{ff}$ is the free-fall time, and $\epsilon_\text{ff} = 0.01$ ([Krumholz2019](https://doi.org/10.1146/annurev-astro-091918-104430)) is the fraction of a cloud's mass that is transformed into stars per cloud free-fall time (also known as star-formation efficiency). 

We have that $t_\text{ff}$ and $\epsilon_\star$ can be written as

$\begin{align}
    t_\text{ff} &= \sqrt{\frac{3\pi}{32 \, G\, \rho_g}} \, , \\
    \epsilon_\star &= s_f(t \rightarrow \infty) = \left. \frac{\rho_s}{\rho} \right\rvert_{t \rightarrow \infty} \, .
\end{align}$

Given that at the beginning of the simulation we don't know $s_f(t \rightarrow \infty)$, we can parametrize it with

$\begin{equation}
    \epsilon_\star = \kappa\,g_f(t_0) \, ,
\end{equation}$
where $\kappa$ is a new dimensionless free parameter that characterizes how much of the gas will be converted to stars, given enough time.

So, we finally have 

$\begin{equation}
    \tau_S = \frac{C_1 \, g_f(t_0)}{\sqrt{\rho_g}} = \frac{C_1 \, g_f(t_0)}{\sqrt{g_f \, \rho(t_0)}} \, ,
\end{equation}$
where we used that $\rho(t) = \rho(t_0)$ thanks to mass conservation, and

$\begin{equation}
    C_1 = \frac{\kappa}{\epsilon_\text{ff}} \, \sqrt{\frac{3\pi}{32 \, G}} \, .
\end{equation}$

We will use $\kappa = 1$, as the simplest choice for the free parameter.
"""

# ╔═╡ 68732d91-805a-4663-9166-f8483213a8d2
begin
	ϵff = 0.01
	κ = 1.0
end;

# ╔═╡ 27281e53-e519-4ad0-af5d-59fb0e208534
C₁ = (κ / ϵff) * sqrt(3π / 32u"G") |> UnitfulAstro.Gyr * UnitfulAstro.Msun^(1/2) * UnitfulAstro.pc^(-3/2)

# ╔═╡ 6503fb74-c34f-40db-afb4-7affd4ceef88
md"""

### Recombination time

From dimensional analysis, we have that the characteristic time for the recombination reaction

$\begin{equation}
    e^- + p^+ \longleftrightarrow H + \gamma \, ,
\end{equation}$

is given by

$\begin{equation}
    \tau_R = \frac{1}{n_e \langle\sigma \, v\rangle} \, , 
\end{equation}$

where $n_e$ is the number fraction of electrons and $\langle\sigma\,v\rangle$ is the recombination rate.

From [Osterbrock2006](http://www.worldcat.org/oclc/60611705) (pg. 22) we have that

$\begin{equation}
    \langle\sigma\,v\rangle \approx \alpha_B(10^4 \, K) = 2.59 \times 10^{-13} \, \mathrm{cm}^3 \, \mathrm{s}^{-1} \, ,
\end{equation}$

where we used the case B recombination, which considers all but one transition channel. A transition directly to the ground state does not result in net recombination, because the emitted photon has enough energy to ionize another atom. Thus it is a better approximation to exclude that case.

The electron density is $n_e = n_i$, where $n_i$ in the number density of ionized Hydrogen atoms, so $n_e$ is essentially the same quantity as $\rho_i(t)$, where the only difference is a conversion factor for the different units, $n_e = n_i = \rho_i(t) / m_p$, where $m_p$ is the proton mass.

We have then

$\begin{equation}
    \tau_R = \frac{C_2}{\rho_i} = \frac{C_2}{i_f \, \rho} \, , 
\end{equation}$
where

$\begin{equation}
	C_2 = \frac{m_p}{\alpha_B(10^4 \, K)} \, .
\end{equation}$
"""

# ╔═╡ f6251e55-f88b-4f53-8449-e30b0bf9ae44
αB = 2.59e-13u"cm^3 * s^-1";

# ╔═╡ 7099a821-a0f0-4931-b6cf-88581e9cff9e
C₂ = u"mp" / αB |> UnitfulAstro.Gyr * UnitfulAstro.Msun * UnitfulAstro.pc^(-3)

# ╔═╡ 4a7eb24b-0874-49a3-9b08-4ffb6a7f0ce7
md"""

### Condensation time

As with the recombination process, we will only consider Hydrogen as a first approximation. There are several channels for the formation of molecular Hydrogen, but the most efficient processes involve the interaction of $H$ atoms on the surface of dust grains, so as before, we have

$\begin{equation}
    \tau_C = \frac{1}{2\,n_\mathrm{dust} \, \langle\sigma v\rangle_\mathrm{dust}}  \, , 
\end{equation}$

where $n_\mathrm{dust}$ is the volume density of dust grains in the ISM, and $\langle\sigma v\rangle_\mathrm{dust}$ is the thermally averaged cross-section for the formation of molecular hydrogen (see [Millán-Irigoyen2020](https://doi.org/10.1093/mnras/staa635) and reference therein).

From [Mollá2017](https://doi.org/10.1093/mnras/stx419) we have that

$\begin{equation}
    n_\mathrm{dust} \, \langle\sigma v\rangle_\mathrm{dust} \approx \frac{Z + Z_\mathrm{eff}}{Z_\odot} \, n_\mathrm{ISM} \, \langle\sigma v\rangle_\odot \, , 
\end{equation}$

where $n_\mathrm{ISM}$ denotes the gas density, $Z$ is the metallicity, $Z_\odot = 0.0134$ ([Asplund2009](https://doi.org/10.1146/annurev.astro.46.060407.145222)) is the solar metallicity, $\langle\sigma v\rangle_\odot = 6 \times 10^{-17} \, \mathrm{cm}^3 \, \mathrm{s}^{-1}$ ([Draine1996](http://doi.org/10.1086/177689)) is the rate of molecular hydrogen formation for $T = 100 \, K$, and $Z_\mathrm{eff} \approx 10^{-3} \, Z_\odot$ ([Glover2007](http://dx.doi.org/10.1086/519445)) is an initial value of metallicity needed to kick start the star formation process, given that the initial abundance of metals and dust grains is zero, and stars only form from molecular clouds. This initial value accounts for all other channels of molecular formation. A detailed study of which would have a minimal impact on the results.

Like before we have $n_\mathrm{ISM} = n_i + n_a + 2\,n_m = \rho_g / m_p$. So, the characteristic time is given by

$\begin{equation}                            \\
    \tau_C = \frac{C_3}{g_f \, \rho}\,\frac{1}{z_f + Z_\mathrm{eff}} \, ,
\end{equation}$

where

$\begin{align}
	g_f &= i_f + a_f + m_f \, , \\
	C_3 &= \frac{Z_\odot \, m_p}{2 \, \langle\sigma v\rangle_\odot} \, .
\end{align}$

Notice that we did the following approximation

$\begin{align}
Z = \frac{\rho_z}{\rho_g} = \frac{z_f \, \rho}{\rho_g} \approx z_f \, .
\end{align}$
"""

# ╔═╡ f2a6676f-457a-476a-9ce7-c336aa9bf47f
begin
	σv = 6e-17u"cm^3 * s^-1"
	Zsun = 0.0134
end;

# ╔═╡ 040e1a8c-97ab-4751-a556-ed936fe58c35
C₃ = Zsun * u"mp" / (2 * σv) |> UnitfulAstro.Gyr * UnitfulAstro.Msun * UnitfulAstro.pc^(-3)

# ╔═╡ 3767c7f9-a0bc-467a-a20a-5e5a266111c7
md"""

## Mass quotients

We define two mass quotients, ionized mass rate per unit of created stellar mass, 

$\begin{equation}
    \eta_\text{ion} = \frac{\dot{M}_\mathrm{ion}}{\psi}  \, , 
\end{equation}$

and disassociated mass rate per unit of created stellar mass, 

$\begin{equation}
    \eta_\text{diss} = \frac{\dot{M}_\mathrm{diss}}{\psi}  \, , 
\end{equation}$

The mass rates $\dot{M}_\mathrm{diss}$ and $\dot{M}_\mathrm{ion}$ can be computed from the photon production rate, for the corresponding energy ranges,

$\begin{align}
    \dot{M}_\mathrm{ion} &= \dot{N}_\mathrm{ion} \, f_\mathrm{ion} \, , \\
    \dot{M}_\mathrm{diss} &= \dot{N}_\mathrm{diss} \, f_\mathrm{diss} \, ,
\end{align}$

where $\dot{N}_\mathrm{ion}$ is the number of ionizing photons produced per unit time (between $0$ and $912\,Å$), $\dot{N}_\mathrm{diss}$ the number of photodissociating photons produced per unit time (in the Lyman–Werner band, $912\,Å$ to $1107\,Å$), and $f_\mathrm{ion}$ and $f_\mathrm{diss}$ are the factors that convert units (proton mass into solar mass) and take into account that the reaction may not be $100\%$ efficient.
"""

# ╔═╡ 127e1dfa-62d8-4721-adc8-cb24c6e9cdcc
md"""
For the ionization reaction, each photon will produce one proton, and we assume $100\%$ efficiency.
"""

# ╔═╡ 005957d6-6f27-4abc-a872-45cf6a032b9f
f_ion = ustrip(1.0u"mp" |> u"Msun")

# ╔═╡ 0fcd2ad5-440c-4128-be21-1f8a354074fe
md"""
For the molecule dissociation reaction, we have to consider the numerical factor given by [Draine1996](https://doi.org/10.1086/177689), where it is shown that dust grains may absorb up to $\sim \! 60$ percent of the photons capable of dissociating hydrogen molecules, and a large fraction of the remaining photons excite different rotational and vibrational states, reducing their dissociation probability to $\sim \! 15$ percent, so we end up with an efficiency factor of $0.4 \times 0.15 = 0.06$.
$f_\mathrm{diss}$ has an extra factor of two with respect to $f_\mathrm{ion}$ because each photon contributes with two protons (from the dissociated molecule) to the atomic gas.
"""

# ╔═╡ f8b02d00-ff30-480e-b5eb-e150e4678c95
f_diss = 0.4 * 0.15 * 2.0 * ustrip(1.0u"mp" |> u"Msun")

# ╔═╡ 44c88ad8-a8c3-45e3-9a56-be3ce5bf66fa
md"""
The number of photons can be computed from $Q(t', Z)$, which is defined as the number of photons produced by a stellar population of one solar mass, of age $t'$ and metallicity $Z$, per unit time. So, we have the relation

$\begin{equation}
    \dot{N}(t) = \int_0^t \psi(t - t') \, Q(t', Z(t - t')) \mathrm{d}t' \, ,
\end{equation}$

where $\psi(t - t')$ is the instantaneous SFR at the moment of birth of the stellar population of current age $t'$, and $Z = Z(t - t')$ is defined as the metallicity at that moment. Because most of the contribution to the integral comes from young blue stars that die in the first $10 \ \mathrm{to} \ 100 \, \mathrm{Myr}$, it is possible to approximate

$\begin{align}
    \psi(t - t') &\approx \psi(t) \, , \\
	Z(t - t') &\approx Z(t) \, .
\end{align}$

So we end up with 

$\begin{equation}
    \eta = f \, \frac{\dot{N}}{\psi} = f \, \int_0^t Q(t', Z) \mathrm{d}t'\, ,
\end{equation}$

The values of $Q$ can be obtained using

$\begin{equation}
    Q(t', Z) = \int_{\lambda_1}^{\lambda_2} \frac{L_\lambda(t', Z)}{E_\lambda} \mathrm{d}\lambda \, ,
\end{equation}$

where $L_\lambda(t', Z)$ is the luminosity per unit of wavelength of a stellar population of one solar mass, of age $t'$ and metallicity $Z$ and $E_\lambda = \lambda / (h \, c)$ is the energy of a photon of wavelength $\lambda$. So integrating 

$\begin{equation}
    \frac{L_\lambda(t', Z)}{E_\lambda} \mathrm{d}\lambda \, ,
\end{equation}$

between the wavelength of interest, we get the number of photons produced per unit time for a stellar population of the given characteristics.

The luminosity will not only depend on the age and metallicity of the population but on the IMF (initial mass function) too, so in principle for each IMF, $t'$ and $Z$ we have a function of luminosity versus $\lambda$.

Using the values from PopStar by [Mollá2009](https://doi.org/10.1111/j.1365-2966.2009.15160.x) we compute a table of $Q_\mathrm{ion}$ and $Q_\mathrm{diss}$ for six IMFs ([Salpeter1955](https://doi.org/10.1086/145971) in two mass ranges, $0.85\,\mathrm{M}_\odot$ to $120\,\mathrm{M}_\odot$ and $0.15\,\mathrm{M}_\odot$ to $100\,\mathrm{M}_\odot$, [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), [Ferrini1990](https://ui.adsabs.harvard.edu/abs/1990A%26A...231..391F) and [Chabrier2003](https://doi.org/10.1086/374879)), six metallicities (0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05) and for ages between $0.1\,\mathrm{Myr}$ and $15\,\mathrm{Gyr}$.
"""

# ╔═╡ 448e1dee-4628-4c14-9d6f-dc165b2e826e
begin
	
	# Raw luminosity data from 
	# https://www.fractal-es.com/PopStar/#download (PopStar2009)
	q_dirs = readdir("./data/luminosity", join = true)
	
	# Regex patterns to extract data from filenames: IMF_mlow_mup_zXXXX_tXXXX
	patterns = [
		r"(?<=spneb_)(.*?)(?=_z)",  # IMF_mlow_mup 
		r".+?(?=_)",                # IMF
		r"(?<=_)(.*?)(?=_)",        # mlow
		r"[^_]+$",                  # mup
		r"(?<=z)(.*?)(?=_)",        # metallicity
		r"(?<=t)(.*)",              # log(age)
	]
	
	# Wavelength range for the photoionization of Hydrogen atoms
	λ_Qion = (0.0u"Å", 912.0u"Å")
	# Wavelength range for the photoionization of Hydrogen molecules
	λ_Qdiss = (912.0u"Å", 1107.0u"Å")
	
	q_data_in_files = DataFrame[]
	
	for dir in q_dirs
		
		files = readdir(dir, join=true)
		
		IMF_mlow_mup = getfield.(match.(patterns[1], basename.(files)), :match)
		IMF = getfield.(match.(patterns[2], IMF_mlow_mup), :match)
		mlow = getfield.(match.(patterns[3], IMF_mlow_mup), :match)
		mup = getfield.(match.(patterns[4], IMF_mlow_mup), :match)
		Zmet = getfield.(match.(patterns[5], basename.(files)), :match)
		ages = getfield.(match.(patterns[6], basename.(files)), :match)
		
		Qion = []
		Qdiss = []
		
		for file in files
			
			data = readdlm(file)
			df = identity.(DataFrame(data, ["λ", "L⋆", "Lneb", "Ltot"]))

			# Wavelength
			df[!, 1] = df[!, 1] .* u"Å"
			# Stellar spectral energy distributions per unit wavelength
			df[!, 2] = df[!, 2] .* 3.82e33u"erg * s^-1 * Å^-1"
			# Nebular spectral energy distributions per unit wavelength
			df[!, 3] = df[!, 3] .* 3.82e33u"erg * s^-1 * Å^-1"
			# Total spectral energy distributions per unit wavelength
			df[!, 4] = df[!, 4] .* 3.82e33u"erg * s^-1 * Å^-1"

			# Spectral energy distribution integration
			let
				λ = @subset(df, λ_Qion[1] .< :λ .< λ_Qion[2])
				integrand = λ[!, 1] .* λ[!, 2] ./ (1u"h" * 1u"c")
				push!(Qion, trapz(λ[!, 1], integrand) |> u"s^-1")
			end
			let
				λ = @subset(df, λ_Qdiss[1] .< :λ .< λ_Qdiss[2])
				integrand = λ[!, 1] .* λ[!, 2] ./ (1u"h" * 1u"c")
				push!(Qdiss, trapz(λ[!, 1], integrand) |> u"s^-1")
			end
			
		end
		
		push!(
			q_data_in_files, 
			identity.(DataFrame(
				"IMF" => uppercase.(IMF),        # Initial mass function
				"mlow" => parse.(Float64, mlow), # Min. mass of the IMF
				"mup" => parse.(Float64, mup),   # Max. mass of the IMF
				"Zmet" => parse.(Float64, "0." .* Zmet),  # Metallicities
				"log(age)" => parse.(Float64, ages),      # Stellar ages
				"Q_ion" => Qion, # Number of ionizing photons per unit time
				"Q_diss" => Qdiss  # Number of dissociating photons per unit time
			))
		)
		
	end
	
	Q_data = sort(vcat(q_data_in_files...), ["IMF", "mlow", "Zmet", "log(age)"])
	
end

# ╔═╡ c7409abf-dc22-429e-ad4d-e2cbd465d454
# Separate the computed Q values by IMF in different variables
begin
	Salpeter1955A = @select(
	    @subset(Q_data, :IMF .== "SAL", :mlow .== 0.85), 
	    $(Not([:IMF, :mlow, :mup])),
    )
	Salpeter1955B = @select(
	    @subset(Q_data, :IMF .== "SAL", :mlow .== 0.15), 
	    $(Not([:IMF, :mlow, :mup])),
    )
	Ferrini1990 = @select(
		@subset(Q_data, :IMF .== "FER"), 
		$(Not([:IMF, :mlow, :mup])),
	) 
    Kroupa2001 = @select(
		@subset(Q_data, :IMF .== "KRO"), 
		$(Not([:IMF, :mlow, :mup])),
	) 
    Chabrier2003 = @select(
		@subset(Q_data, :IMF .== "CHA"), 
		$(Not([:IMF, :mlow, :mup])),
	)  
    Q_imfs = [
		Salpeter1955A, 
		Salpeter1955B,
		Ferrini1990,
		Kroupa2001, 
		Chabrier2003,
	]
end;

# ╔═╡ aa5e9990-db35-4a91-912e-f839daf6c686
begin
	# Names of the columns for some dataframes
	col_names = [
		"Z", 
		"Kroupa2001", 
		"Ferrini1990", 
		"Salpeter1955A", 
		"Salpeter1955B",
		"Chabrier2003",
	]

	# Metallicities we have data for
    Q_metals = unique(Kroupa2001[!, "Zmet"])

	# Maximum stellar age used to compute ηion and ηdiss: 16 Gyr
	MAX_AGE = log10(16e9) 

	# Metallicity used to compute ηion and ηdiss: 0.008
	METAL = 0.008
end;

# ╔═╡ 342c1ad8-338e-44e2-adec-7638fe1767a2
# Save the η values in text files
# Each file is named after the metallicity and is 
# stored within the folders ./data/eta_data/`IMF`/
for (name, imf) in zip(col_names[2:end], Q_imfs)
	for z in Q_metals
		dir = mkpath("./data/eta_data/$name/")
			
		eta_data = Array{Float64}(undef, 0, 3)
		q_df = @select(@subset(imf, :Zmet .== z), $(Not(:Zmet)))

		for index_age in 5:size(q_df, 1)
				
			ages = exp10.(q_df[1:index_age,  "log(age)"]) * u"yr"
		    QH = q_df[1:index_age, "Q_ion"]
			QHII = q_df[1:index_age, "Q_diss"]
				
			eta_data = vcat(
				eta_data, 
				[
					q_df[index_age,  "log(age)"];;
					ustrip(Unitful.NoUnits, trapz(ages, QH) * f_ion);;
					ustrip(Unitful.NoUnits, trapz(ages, QHII) * f_diss)
				],
			)
		end
			
		writedlm(joinpath(dir, string(z)), eta_data)
	end
end

# ╔═╡ 7788b98a-5bec-4b6d-82d9-2c272e2255a7
md"""### Ionization efficiency of Hydrogen atoms ($\eta_\mathrm{ion}$)"""

# ╔═╡ 7af2d4a6-a304-404d-8cbe-eeddb80beba6
md"""
Now that we have all the $Q_\mathrm{ion}$ values we will integrate them (up to the maximum  stellar age we have data for ~$15 \, \mathrm{Gyr}$) for each IMF and each value of $Z$, resulting in the following table
"""

# ╔═╡ aeb72f0e-2252-486a-b79b-9d8cc6e5f962
begin
    table_ion = Array{Float64}(undef, 6, 0)
	
    for imf in Q_imfs
        sol = Float64[]
        for met in Q_metals
			sub_df = @subset(imf, :Zmet .== met, $("log(age)") .< MAX_AGE)
			ages = exp10.(sub_df[!, "log(age)"]) * u"yr"
    	    QH = sub_df[!, "Q_ion"]
    	    append!(sol, uconvert(Unitful.NoUnits, trapz(ages, QH) * f_ion))
        end
        global table_ion = hcat(table_ion, sol)
    end
	
    table_ion = hcat(Q_metals, table_ion)
	eta_ion = DataFrame(table_ion, col_names)
end

# ╔═╡ e83337bd-8c2d-4a9a-bd8b-7f8201cf67ad
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"Z", 
		ylabel = L"\eta_\mathrm{ion}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		xscale = log10,
	)

	for col in col_names[2:end]
		scatterlines!(
			ax, 
			Q_metals, eta_ion[!, col], 
			linewidth = 3, 
			label = col,
		)
	end

	axislegend(position = :rt)

	f
end

# ╔═╡ b8beaaa6-b018-4ca8-b19d-a918dc761707
md"""#### Ionization efficiency as a function of the maximum stellar age"""

# ╔═╡ 43eb3af9-86c5-49e9-af0e-3270e3df493e
md"""
At fix metallicity (Z = 0.008).
"""

# ╔═╡ 05d22c2d-f76f-4931-8e21-6d31e9ab178e
begin
	AGES = [
		log10(1e6), log10(2e6), log10(3e6), log10(4e6), log10(5e6), log10(6e6), log10(7e6), log10(8e6), log10(9e6), log10(10e6), log10(100e6), log10(1e9), log10(10e9)
	]
    table_z_ion = Array{Float64}(undef, length(AGES), 0)
	
    for imf in Q_imfs
        sol = Float64[]
        for age in AGES
			sub_df = @subset(imf, :Zmet .== METAL, $("log(age)") .< age)
			ages = exp10.(sub_df[!, "log(age)"]) * u"yr"
    	    QH = sub_df[!, "Q_ion"]
    	    append!(sol, uconvert(Unitful.NoUnits, trapz(ages, QH) * f_ion))
        end
    	global table_z_ion  = hcat(table_z_ion, sol)
	end

	table_z_ion = hcat(AGES, table_z_ion)
	eta_z_ion = DataFrame(table_z_ion, ["log(age)", col_names[2:end]...])
end

# ╔═╡ ea0bada1-4359-4c2b-9dc8-d91b6ebd5686
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"\mathrm{Stellar} \,\, \mathrm{age}", 
		ylabel = L"\eta_\mathrm{ion}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		titlesize = 25,
		xscale = log10,
		title = L"Z = 0.008",
	)

	for col in col_names[2:end]
		scatterlines!(
			ax, 
			AGES, eta_z_ion[!, col], 
			linewidth = 3, 
			label = col,
		)
	end

	axislegend(position = :rt)

	f
end

# ╔═╡ 00e1f4f0-6a6a-48dc-81d5-99f356aa3410
md"""
At fix IMF (Kroupa et al. 2001).
"""

# ╔═╡ 0f212d88-6b3a-4f96-8f05-d5b0f9943fe3
begin
    table_imf_ion = Array{Float64}(undef, length(AGES), 0)
	
    for met in Q_metals
        sol = Float64[]
        for age in AGES
			sub_df = @subset(Q_imfs[1], :Zmet .== met, $("log(age)") .< age)
			ages = exp10.(sub_df[!, "log(age)"]) * u"yr"
    	    QH = sub_df[!, "Q_ion"]
    	    append!(sol, uconvert(Unitful.NoUnits, trapz(ages, QH) * f_ion))
        end
    	global table_imf_ion  = hcat(table_imf_ion, sol)
	end

	table_imf_ion = hcat(AGES, table_imf_ion)
	eta_imf_ion = DataFrame(table_imf_ion, ["log(age)", string.(Q_metals)...])
end

# ╔═╡ 3c26c7f7-0d0c-4fe3-a071-ab76c7904659
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"\mathrm{Stellar} \,\, \mathrm{age}", 
		ylabel = L"\eta_\mathrm{ion}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		titlesize = 25,
		xscale = log10,
		title = L"Kroupa et al. $2001$"
	)

	for met in Q_metals
		scatterlines!(
			ax, 
			AGES, eta_imf_ion[!, string(met)], 
			linewidth = 3, 
			label = string(met),
		)
	end

	axislegend(position = :lt)

	f
end

# ╔═╡ fa38bb0f-c807-48b2-8bf2-33e457b53435
md"""### Approximations of $\eta_\mathrm{ion}$"""

# ╔═╡ e879e854-a9c3-4762-9963-dfa8960f2dc5
md"""#### Fitting"""

# ╔═╡ 66f233f2-76d8-45be-9aa2-831ff7269d96
md"""
Given that the result of $\eta_\mathrm{ion}$ for a given IMF is a function of Z, we can fit a quadratic function to the log-log data

$\begin{equation}
    \log_{10}(\eta) = p_1 \log_{10}(Z)^2 + p_2 \log_{10}(Z) + p_3 \, ,
\end{equation}$
where $p_1$, $p_2$ and $p_3$ are the parameters to be fitted.
"""

# ╔═╡ d7069cdb-9a99-4194-93ea-03b5b82ddc89
# Fitting model
begin
	@. model(x, p) = exp10(p[1] * log10(x)^2 + p[2] * log10(x) + p[3])
	
	xdata = Q_metals
	p0 = [-1.0, -1.0, 1.0]
end;

# ╔═╡ 4525cc1d-5338-400e-840a-47c069014ec6
md"""The fitted parameters for each IMF are"""

# ╔═╡ 0e350f3a-534e-4b7b-b473-9f3ebefb399e
begin
	ion_coeff = DataFrame("coeff" => ["p₁", "p₂", "p₃"])
	
	for (ydata, colname) in zip(eachcol(eta_ion[!, 2:end]), col_names[2:end])
		fit = curve_fit(model, xdata, ydata, p0)
		ion_coeff[!, colname] = coef(fit)
	end

	ion_coeff
end

# ╔═╡ fd2f6a01-9125-49fb-8615-b2e5deccd107
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"Z", 
		ylabel = L"\eta_\mathrm{ion}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		titlesize = 25,
		yscale = log10,
		xscale = log10,
		title = "Fitting"
	)

	for (i, col) in enumerate(col_names[2:end])
		scatter!(
			ax, 
			Q_metals, eta_ion[!, col], 
			linewidth = 3, 
			label = col,
		)
		lines!(
			ax, 
			Q_metals,  model(Q_metals, ion_coeff[!, i + 1]), 
			linewidth = 3,
			color = :red,
		)
	end	

	axislegend(position = :lb)

	f
end

# ╔═╡ 930fe232-f4ff-4331-9752-4e4ec0e66009
md"""#### Interpolation"""

# ╔═╡ 386334a4-60bf-47de-bf72-516f042b1407
md"""Alternatively we can interpolate linearly between the known values, setting a flat constant (equal to the las know value) for metallicities outside the range."""

# ╔═╡ 9b0f3856-bee3-485e-9049-1a082c86e571
# Interpolation
itp_ion = Dict(
	colname => LinearInterpolation(xdata, ydata, extrapolation_bc = Flat()) for 
	(ydata, colname) in zip(eachcol(eta_ion[!, 2:end]), col_names[2:end])
);

# ╔═╡ d5ba04de-e2e8-44af-9b60-6a47b782248e
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"Z", 
		ylabel = L"\eta_\mathrm{ion}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		titlesize = 25,
		yscale = log10,
		xscale = log10,
		title = "Interpolation"
	)

	z_vals = 10 .^ range(-5, stop=-1, length=100)

	for col in col_names[2:end]
		scatter!(
			ax, 
			Q_metals, eta_ion[!, col], 
			linewidth = 3, 
			label = col,
		)
		lines!(
			ax, 
			z_vals,  itp_ion[col](z_vals), 
			linewidth = 3,
			color = :red,
		)
	end	

	axislegend(position = :lb)

	f
end

# ╔═╡ 396aa4fe-5d65-4852-9cb2-03c654201f6e
md"""### Photodissociation efficiency of $\mathrm{H}_2$ molecules ($\eta_\mathrm{diss}$)"""

# ╔═╡ 399a13c4-235b-40e7-8a2f-5affd889c014
md"""
As we did for $\eta_\mathrm{ion}$, we integrate the $Q_\mathrm{diss}$ values (up to the maximum  stellar age we have data for ~$15 \, \mathrm{Gyr}$) for each IMF and each value of $Z$, resulting in the following table
"""

# ╔═╡ f6c88e38-5f80-4ab7-afef-4f3249af8723
begin
    table_diss = Array{Float64}(undef, 6, 0)
	
    for imf in Q_imfs
        sol = Float64[]
        for met in Q_metals
			sub_df = @subset(imf, :Zmet .== met, $("log(age)") .< MAX_AGE)
	    	ages = exp10.(sub_df[!, "log(age)"]) * u"yr"
    	    QHII = sub_df[!, "Q_diss"]
    	    append!(sol, uconvert(Unitful.NoUnits, trapz(ages, QHII) * f_diss))
        end
        global table_diss = hcat(table_diss, sol)
    end
	
    table_diss = hcat(Q_metals, table_diss)
	eta_diss = DataFrame(table_diss, col_names)
end

# ╔═╡ ff8a6018-dd86-4b09-a122-a72e0dfa7013
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"Z", 
		ylabel = L"\eta_\mathrm{diss}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		yscale = log10,
		xscale = log10,
	)

	for col in col_names[2:end]
		scatterlines!(
			ax, 
			Q_metals, eta_diss[!, col], 
			linewidth = 3, 
			label = col,
		)
	end

	axislegend(position = :lb)

	f
end

# ╔═╡ 2c58b32c-e731-467b-b051-7063b3d3e341
md"""#### Photodissociation efficiency as a function of the maximum stellar age"""

# ╔═╡ 2b7d684e-3db3-4966-9993-0446a7db7edb
md"""
At fix metallicity (Z = 0.008).
"""

# ╔═╡ 0a465f30-6904-4195-830e-21cfb00fb63a
begin
    table_z_diss = Array{Float64}(undef, length(AGES), 0)
	
    for imf in Q_imfs
        sol = Float64[]
        for age in AGES
			sub_df = @subset(imf, :Zmet .== METAL, $("log(age)") .< age)
			ages = exp10.(sub_df[!, "log(age)"]) * u"yr"
    	    QHII = sub_df[!, "Q_diss"]
    	    append!(sol, uconvert(Unitful.NoUnits, trapz(ages, QHII) * f_diss))
        end
    	global table_z_diss  = hcat(table_z_diss, sol)
	end

	table_z_diss = hcat(AGES, table_z_diss)
	eta_z_diss = DataFrame(table_z_diss, ["log(age)", col_names[2:end]...])
end

# ╔═╡ d7871e31-d841-4951-aec5-d6f2915d0ceb
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"\mathrm{Stellar} \,\, \mathrm{age}", 
		ylabel = L"\eta_\mathrm{diss}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		titlesize = 25,
		xscale = log10,
		title = L"Z = 0.008",
	)

	for col in col_names[2:end]
		scatterlines!(
			ax, 
			AGES, eta_z_diss[!, col], 
			linewidth = 3, 
			label = col,
		)
	end

	axislegend(position = :lt)

	f
end

# ╔═╡ 9b738a7d-1a85-460b-9585-6ee4f32b41ae
md"""
At fix IMF ([Kroupa et al. 2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x)).
"""

# ╔═╡ afdb1375-ac94-4dcc-8fcb-c732e6029b83
begin
    table_imf_diss = Array{Float64}(undef, length(AGES), 0)
	
    for met in Q_metals
        sol = Float64[]
        for age in AGES
			sub_df = @subset(Q_imfs[1], :Zmet .== met, $("log(age)") .< age)
			ages = exp10.(sub_df[!, "log(age)"]) * u"yr"
    	    QHII = sub_df[!, "Q_diss"]
    	    append!(sol, uconvert(Unitful.NoUnits, trapz(ages, QHII) * f_diss))
        end
    	global table_imf_diss  = hcat(table_imf_diss, sol)
	end

	table_imf_diss = hcat(AGES, table_imf_diss)
	eta_imf_diss = DataFrame(table_imf_diss, ["log(age)", string.(Q_metals)...])
end

# ╔═╡ 1e2f1eff-3295-4eae-ac83-569e9f8211fd
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"\mathrm{Stellar} \,\, \mathrm{age}", 
		ylabel = L"\eta_\mathrm{diss}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		titlesize = 25,
		xscale = log10,
		title = L"Kroupa et al. $2001$"
	)

	for met in Q_metals
		scatterlines!(
			ax, 
			AGES, eta_imf_diss[!, string(met)], 
			linewidth = 3, 
			label = string(met),
		)
	end

	axislegend(position = :lt)

	f
end

# ╔═╡ 67dfb1b3-4818-445f-bfe2-16d442506567
md"""### Approximations of $\eta_\mathrm{diss}$"""

# ╔═╡ 3da12a3d-bc26-4115-b769-78841499e434
md"""#### Fitting"""

# ╔═╡ 64fda1f0-8fcf-4898-b3dd-883f151f69b3
md"""Using the same fitting function as for $\eta_\mathrm{ion}$, we get"""

# ╔═╡ 76136e33-14f7-4c5e-84e9-970372c6c01a
begin
	diss_coeff = DataFrame("coeff" => ["p1", "p2", "p3"])
	
	for (ydata, colname) in zip(eachcol(eta_diss[!, 2:end]), col_names[2:end])
		fit = curve_fit(model, xdata, ydata, p0)
		diss_coeff[!, colname] = coef(fit)
	end

	diss_coeff
end

# ╔═╡ 1fce7934-0ce4-4059-ba8e-691f1c505f4d
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"Z", 
		ylabel = L"\eta_\mathrm{diss}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		titlesize = 25,
		yscale = log10,
		xscale = log10,
		title = "Fitting"
	)

	for (i, col) in enumerate(col_names[2:end])
		scatter!(
			ax, 
			Q_metals, eta_diss[!, col], 
			linewidth = 3, 
			label = col,
		)
		lines!(
			ax, 
			Q_metals,  model(Q_metals, diss_coeff[!, i + 1]), 
			linewidth = 3,
			color = :red,
		)
	end	

	axislegend(position = :lb)

	f
end

# ╔═╡ e9eca322-1466-4d18-af62-28d3edd42cfc
md"""#### Interpolation"""

# ╔═╡ 68502224-d186-464b-8d35-1ce9ad0a9994
# Interpolation
itp_diss = Dict(
	colname => LinearInterpolation(xdata, ydata, extrapolation_bc = Flat()) for 
	(ydata, colname) in zip(eachcol(eta_diss[!, 2:end]), col_names[2:end])
);

# ╔═╡ f017b433-05dd-48bc-b0b4-f70a39100b2d
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"Z", 
		ylabel = L"\eta_\mathrm{diss}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		titlesize = 25,
		yscale = log10,
		xscale = log10,
		title = "Interpolation"
	)

	z_vals = 10 .^ range(-5, stop=-1, length=100)

	for col in col_names[2:end]
		scatter!(
			ax, 
			Q_metals, eta_diss[!, col], 
			linewidth = 3, 
			label = col,
		)
		lines!(
			ax, 
			z_vals,  itp_diss[col](z_vals), 
			linewidth = 3,
			color = :red,
		)
	end	

	axislegend(position = :lb)

	f
end

# ╔═╡ 533b3cd0-c1f6-4ecd-b196-4ed35bf77135
md"""

## Mass recycling

There are two mass recycling parameters, $R$ which is defined as the mass fraction of a stellar population that is returned to the ISM under the instantaneous recycling hypothesis (stars under certain mass live forever, and stars above that mass die instantly), and $Z_\mathrm{SN}$ which is the fraction of the returned gas that is composed of metals (the rest is assumed to be ionized gas). 
A stellar yield model gives the amount (as a fraction of the stellar mass) of each modeled element that is returned to the ISM by stars with masses between $m$ and $m + \mathrm{d}m$, so they can be used to compute

$\begin{equation}
	R = \dfrac{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \mathrm{d}m}{\int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \mathrm{d}m} \, ,
\end{equation}$

where $\phi(m)$ is the initial mass function (ISM), $m_\mathrm{low}$ and $m_\mathrm{high}$ are the extremes in the mass range of the ISM used, $m_\mathrm{ir}$ is the mass limit for the instantaneous mass recycling hypothesis, and $m_\mathrm{rem}(m)$ is the remnant stellar mass given by the yield model.

Using the same notation we can calculate $Z_\mathrm{SN}$ as

$\begin{equation}
	Z_\mathrm{SN} = \dfrac{\int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, f_Z \, \phi(m) \mathrm{d}m}{\int_{m_\mathrm{ir}}^{m_\mathrm{high}} (m - m_\mathrm{rem}(m)) \, \phi(m) \mathrm{d}m} \, ,
\end{equation}$

where $f_Z$ is the mass fraction of the stellar mass that is returned to the ISM as metals.

Notice that the denominator in the expression for $R$ is the total mass of the stellar population modeled by $\phi(m)$, so is only action as a normalization given that the IMF is always defined except for a global constant.

Some traditional choices for the masses are $m_\mathrm{ir} = 8 \, M_\odot$, $m_\mathrm{low} = 0.08 \, M_\odot$ (the limit for hydrogen fusion), and $m_\mathrm{low} = 100 \, M_\odot$ (an order of magnitude for the upper limit of validity of the yield models, given by experimental limitations).

[Ascasibar2015](https://doi.org/10.1093/mnras/stv098) got $R \approx 0.18$ and $Z_{SN} \approx 0.09$, using the yield model of [Woosley1995](https://doi.org/10.2172/115557) and the IMF of [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x), even though it is not clear which mass limits were used.
"""

# ╔═╡ 97e7f37b-9494-4ae9-a076-77b90a974a81
begin
	model_names = Dict(
		"Woosley1995" => "WOW",
		"Portinari1998" => "PCB", 
		"Chieff2004" => "CLI", 
		"Kobayashi2006" => "KOB",
		"Heger2010" => "HEG",
		"Limongi2012" => "LIM", 
	)
	md"""
	We will consider the stellar yields models by
	
	 -  "WOW" --> [Woosley et al. 1995](https://doi.org/10.2172/115557)
	 -  "PCB" --> [Portinari et al. 1998](https://ui.adsabs.harvard.edu/abs/1998A&A...334..505P)
	 -  "CLI" --> [Chieff et al. 2004](https://doi.org/10.1086/392523)
	 -  "KOB" --> [Kobayashi et al. 2006](https://doi.org/10.1086/508914)
	 -  "HEG" --> [Heger et al. 2010](https://doi.org/10.1088/0004-637X/724/1/341)
	 -  "LIM" --> [Limongi et al. 2012](https://doi.org/10.1088/0067-0049/199/2/38)
	
	compiled by [Mollá2015](https://doi.org/10.1093/mnras/stv1102), which are summarized in the following table, where
	  - `model`: Stellar yield model
	  - `s_Z`: Metallicity of the stellar population modeled by the IMF
	  - `s_m`: Stellar mass
	  - `m_rem`: Remnant mass, after stellar death
	  - `zf_rem`: Fraction of the stellar mass ejected as metals to the ISM
	"""
end

# ╔═╡ be85ba3b-5439-4cf3-bb14-d24d61a283c3
begin
	# Raw stellar yields from 
	# Mollá et al. 2015 (https://doi.org/10.1093/mnras/stv1102)
	sy_files = readdir("./data/stellar_yields", join = true)
	
	sy_data = DataFrame[]
	
	for file in sy_files
		
		data = readdlm(file, skipstart = 1)
		df = identity.(
			DataFrame(
				data, 
				["s_Z", "s_m", "H", "D", "He3", "He4", "C12", "O16", "Ne20", "Mg24", "Si28", "S32", "Ca40", "Fe56", "m_rem", "C13s", "N14s"],
			),
		)

		name = uppercase(getindex(getproperty(
			match(r"^(.+?)(\.[^.]*$|$)", basename(file)), 
			:captures,
		), 1))
		
		# Stellar mass
		df[!, 2] = df[!, 2] .* u"Msun"
		# Remnant mass
		df[!, 15] = df[!, 15] .* u"Msun"

		insertcols!(df, 1, "model" => fill(name, length(df[!, 1])))

		@transform! df begin
           :zf_rem = (
			   :C12 .+ :O16 .+ :Ne20 .+ :Mg24 .+ :Si28 .+ :S32 
			   .+ :Ca40 .+ :Fe56 .+ :C13s .+ :N14s .+ (1 .- :m_rem ./ :s_m) .* :s_Z
		   )
       end

		select!(
			df, 
			Not(["H", "D", "He3", "He4", "C12", "O16", "Ne20", "Mg24", 
				"Si28", "S32", "Ca40", "Fe56", "C13s", "N14s"]),
		)

		push!(sy_data, df)
	end
	
	sy_data = sort(vcat(sy_data...), ["model", "s_Z", "s_m"])
	
	# Metallicities
	sy_metallicities = unique(sy_data[!, "s_Z"])
	# Masses
	sy_masses = unique(sy_data[!, "s_m"]) 

	# Columns:
	#   model: Stellar yield model
	#   s_Z: Metallicity of the stellar population modeled by the IMF
	#   s_m: Stellar mass
	#   m_rem: Remnant mass, after stellar death
	#   zf_rem: Fraction of the stellar mass ejected as metals to the ISM
	sy_data
end

# ╔═╡ d1e89b59-bc6f-46fd-a7cd-126fad530916
md"""
With the choice $m_\mathrm{low} = 0.08 \, M_\odot$, $m_\mathrm{high} = 100 \, M_\odot$, and $m_\mathrm{ir} = 8 \, M_\odot$ we get the following $R$ and $Z_\mathrm{SN}$,
"""

# ╔═╡ 49d39360-3609-407c-bfee-c46e7485727a
md"""
Averaging over stellar metallicity (`s_Z` column), and varying $m_\mathrm{high}$ between $30 \, M_\odot$ and $100 \, M_\odot$, we get the following plots, where the error bars indicate the standard deviation for varying the stellar metallicity.
We choose for the plots the [Kroupa2001](https://doi.org/10.1046/j.1365-8711.2001.04022.x) IMF and the [Woosley et al. 1995](https://doi.org/10.2172/115557) yield model.
"""

# ╔═╡ 9666bdc8-cbc0-4757-9bd8-a76477c252eb
md"# Implementation"

# ╔═╡ ca9a233b-d3ca-4a76-a3d8-f29884ac9484
md"## Constants"

# ╔═╡ e2e4ae4f-dcdc-4999-88f2-853378be859a
md"## Equations"

# ╔═╡ 2be98f2a-57a2-4f53-ad53-b1c0e9e9aafa
md"""
Each equation has units of $\mathrm{Gyr}^{-1}$, the parameter $\rho_0$ has to have units of $\mathrm{M_\odot} \, \mathrm{pc}^{-3}$, and the parameter $g_0$ should be dimensionless.

```
Ionized gas fraction:       i(t) / ρ -> y[1]
Atomic gas fraction:        a(t) / ρ -> y[2]
Molecular gas fraction:     m(t) / ρ -> y[3]
Metal fraction:             z(t) / ρ -> y[4]
Stellar fraction:           s(t) / ρ -> y[5]	

where ρ = i(t) + a(t) + m(t) + s(t)
```   

"""

# ╔═╡ 01c985a1-4d3d-4799-8e9e-5e25390a9fbb
md"""## Jacobian"""

# ╔═╡ ceec9b81-6a43-4bb9-bb74-b309ef4c3037
md"# Functions"

# ╔═╡ 4607856c-7472-4131-a2ee-29f7150f5cb4
md"## Integration"

# ╔═╡ 4cfe1c80-c67e-4dd3-825b-d893800d68c0
md"## Density PDF"

# ╔═╡ d7ba9e0c-5cfa-4176-adff-12cb8e20679b
md"### Parameters for the density PDF"

# ╔═╡ 82e78dc9-b89e-48d9-9f70-6f3238dfd196
Base.@kwdef struct Params
	# Power law slope
	α::Float64
	# Dimensionless turbulent forcing parameter
	b::Float64
	# Mach number
	Ms::Float64
	# (min, max) values of s = log(ρ/ρ₀)
	deviation::NTuple{2,Float64}
	# Number of subdivisions of the variable (s = log(ρ/ρ₀) or f = ρ/ρ₀)
	divisions::Int64
	# Total initial mean density in Mₒ pc^(-3)
	ρ₀::Float64        
end;

# ╔═╡ 7a2987ef-d37e-4c7a-aaa8-8186694bea88
function mass_fraction(
	ρ_PDF::Union{Nothing,Function}, 
	params::Params,
	log_var::Bool,
)::NTuple{2, Vector{Float64}}

	if params.divisions == 1
		return [1,], [log_var ? 0 : 1.0,]
	end

	# Which variable will be used
	# log_var == true:  s = log(ρ/ρ₀)
	# log_var == false: f = ρ/ρ₀
	dev = log_var ? params.deviation : exp.(params.deviation)
			
	# Step in the range of the variable (s = log(ρ/ρ₀) or f = ρ/ρ₀)
    step = (dev[2] - dev[1]) / params.divisions

	# Values of the variable (s = log(ρ/ρ₀) or f = ρ/ρ₀)
	points = [dev[1] + step * (i - 0.5) for i in 1:params.divisions]

	# Fractions of mass within each division
	mass_f = [quadgk(
		x -> ρ_PDF(x, params), 
		log_var ? point - (step / 2) : log(point - (step / 2)), 
		log_var ? point + (step / 2) : log(point + (step / 2)), 
		order = 10,
	)[1] for point in points]

	return mass_f, points
	
end;

# ╔═╡ 34c053bf-0b4f-45c4-bb79-7e5e89a26060
md"### Density PDF by Burkhart et al. (2018)"

# ╔═╡ 768f8ffb-a08b-4498-97e4-1a3a866e69c7
function pBurkhart2018(s, params)
	
	b = params.b
	Ms = params.Ms
	α = params.α
	
	σs2 = log(1 + b^2 * Ms^2)
	s0 = -0.5 * σs2
	st = (α - 0.5) * σs2
	C = exp((α - 1) * 0.5 * α * σs2) / sqrt(2π * σs2)
	N = 1 / ((C * exp(-α * st)) / α + 0.5 + 0.5 * erf((2 * st + σs2) / sqrt(8 * σs2)))
	
	if s < st
		return (N / sqrt(2π * σs2)) * exp(-((s - s0)^2) / (2 * σs2))
	else
		return N * C * exp(-α * s)
	end
	
end;

# ╔═╡ 8864b4d3-6a9f-4e7b-8cd6-ed32a0116f4a
md"### Density PDF by Krumholz et al. (2005)"

# ╔═╡ aad5e227-67a5-49a0-a79d-24160d3ebe06
function pKrumholz2005(s, params)
	
	b = params.b
	Ms = params.Ms
	
	σs2 = log(1 + b^2 * Ms^2)
	s0 = -0.5 * σs2
	
	return exp(-((s - s0)^2) / (2 * σs2)) / sqrt(2π * σs2)
	
end;

# ╔═╡ d6502bbe-ea58-471d-a748-d9c8f37ce61e
# Save the integrated version of the density PDF function 
# in the file ./data/rho_pdf.txt
begin
	const DIVISIONS = 10
	const DEVIATION = [-4, 6]
	const FUNCTION = pBurkhart2018
	const PARAMS = Params(2, 0.5, 10.0, (DEVIATION...,), DIVISIONS, 0)

	step = (DEVIATION[2] - DEVIATION[1]) / DIVISIONS
	first_s = DEVIATION[1] + (step / 2)

	open("./data/rho_pdf.txt", "w") do file
		write(file, "#define DIVISIONS $DIVISIONS\n")
		write(file, "/* Probability density function of the interstellar gas density */\n")
		write(file, "static const double PDF[] = {\n")
		for i in 1:DIVISIONS
			s = first_s + step * (i - 1)
			int_val = quadgk(
				x -> FUNCTION(x, PARAMS), 
				s - (step / 2), 
				s + (step / 2), 
				order = 10,
			)[1]
			write(file, "\t$int_val,\n")
		end
		write(file, "};\n")

		write(file, "/* Density factor: rho = rho0 * F_RHO */\n")
		write(file, "static const double F_RHO[] = {\n")
		for i in 1:DIVISIONS
			f = exp(first_s + step * (i - 1))
			write(file, "\t$f,\n")
		end
		write(file, "};")
	end
end;

# ╔═╡ b3a260b6-eb31-43a0-9fd6-60a507984319
md"""
## Initial mass functions (IMF)

The initial mass function $\phi(m)$ gives the number of stars between mass $m$ and $m + \mathrm{d}m$, for a given population of total mass $M$, given by the relation

$\begin{equation}            
    M = \int_{m_\mathrm{low}}^{m_\mathrm{high}} m \, \phi(m) \, \mathrm{d}m \, ,
\end{equation}$

which allows to normalize $\phi(m)$ for the population of total mass $M$, within the range $[m_\mathrm{low}, m_\mathrm{high}]$.

There are many models for $\phi(m)$, but one of the most common is a simple power law

$\begin{equation}            
    \phi(m) = A \, m^{-\alpha}\, ,
\end{equation}$

where it is generally assumed that $[m] = M_\odot$ and $A$ is the normalization constant.

The following implementations don't have a specific choice for normalization (when posible $A = 1$), so they have to be multiplied by a constant if one wants a given value of $M$.
"""

# ╔═╡ 5ba3a0c1-6107-45a1-9b1d-5c323b9a7145
begin
	# Salpeter 1955 (https://doi.org/10.1086/145971)
	# This model is valid for 0.4 <= m / M⊙ <= 10
	ϕSAL(m::Float64)::Float64 = m^(-2.35)
	
	# Miller et al. 1979 (https://doi.org/10.1086/190629)
	# This model is valid for 0.1 <= m / M⊙ <= 62
	const C1_MIL = 1.09
	const C2_MIL = -1.02
	ϕMIL(m::Float64)::Float64 = m^(-1) * exp(-(log10(m) - C2_MIL)^2 / (1 / C1_MIL))
	
	# Ferrini et al. 1990 (https://ui.adsabs.harvard.edu/abs/1990A%26A...231..391F)
	# Ferrini et al. 1992 (https://doi.org/10.1086/171066)
	# From the papers it is not clear the range of validity for this model, but it is 
	# generaly accepted that no model is valid outside 0.072 <= m / M⊙ <= 100
	ϕFER(m::Real)::Real = m^(-0.52) * exp10(
		-sqrt(0.73 + log10(m) * (1.92 + 2.07 * log10(m)))
	)
	
	# Kroupa 1993 (https://doi.org/10.1093/mnras/262.3.545)	
	# This model is valid for m / M⊙ >= 0.072
	function ϕKRO_93(m::Real)::Real
		if m < 0.5
			return m^(-1.2)
		elseif 0.5 <= m < 1
			return 0.5 * (m^-2.2)
		else
			return 0.5 * (m^-2.7)
		end
	end
	
	# Kroupa 2001 (https://doi.org/10.1046/j.1365-8711.2001.04022.x)
	# This model is valid for m / M⊙ >= 0.072
	function ϕKRO_01(m::Real)::Real
		if m < 0.08
			return m^-0.3
		elseif 0.08 <= m < 0.5
			return 0.08 * (m^-1.3)
		else
			return 0.0386 * (m^-2.35)
		end
	end
	
	# Chabrier 2003 (https://doi.org/10.1086/374879)
	# This model is valid for m / M⊙ <= 10 
	# (above m = 1 M⊙ uses Salpeter (1955) results)
	function ϕCHA(m::Real)::Real
		if m <= 1
			return m^(-1) * exp(-(log10(m) - log10(0.22))^2 / (2 * 0.57^2))
		else
			return 0.514 * m^(-2.3)
		end		
	end	
	
	# Weidner 2005 (https://doi.org/10.1086/429867)
	# This model is valid for m / M⊙ >= 0.072
	function ϕWEI(m::Real)::Real
		if m < 0.08
			return m^(-0.3)
		elseif 0.08 <= m < 0.5
			return 0.08 * m^(-1.3)
		elseif 0.5 <= m < 1
			return 0.0386 * m^(-2.35)
		else
			return 0.0386 * m^(-2.7)
		end
	end

	# Millán-Irigoyen et al. 2020 (https://doi.org/10.1093/mnras/staa635)
	# This model is valid for 0.1 <= m / M⊙ <= 50
	function ϕMILLA(m::Real)::Real
		if m < 0.5
			return m^(-1.3)
		else
			return 0.5 * m^(-2.3)
		end
	end

	ϕSAL(m::Quantity)::Float64 = ϕSAL(ustrip(u"Msun", m))
	ϕMIL(m::Quantity)::Float64 = ϕMIL(ustrip(u"Msun", m))
	ϕFER(m::Quantity)::Float64 = ϕFER(ustrip(u"Msun", m))
	ϕKRO_93(m::Quantity)::Float64 = ϕKRO_93(ustrip(u"Msun", m))
	ϕKRO_01(m::Quantity)::Float64 = ϕKRO_01(ustrip(u"Msun", m))
	ϕCHA(m::Quantity)::Float64 = ϕCHA(ustrip(u"Msun", m))
	ϕWEI(m::Quantity)::Float64 = ϕWEI(ustrip(u"Msun", m))
	ϕMILLA(m::Quantity)::Float64 = ϕMILLA(ustrip(u"Msun", m))

	imf_funcs = Dict(
		"Salpeter1955" => ["SAL", ϕSAL],
		"Miller1979" => ["MIL", ϕMIL],
		"Ferrini1990" => ["FER", ϕFER],
		"Kroupa1993" => ["KRO93", ϕKRO_93],
		"Kroupa2001" => ["KRO01", ϕKRO_01],
		"Chabrier2003" => ["CHA", ϕCHA],
		"Weidner2005" => ["WEI", ϕWEI],
		"Millán-Irigoyen2020" => ["MILLA", ϕMILLA],
	)
end;

# ╔═╡ 946d007d-abd7-4cd3-9789-e77b1ad6ebf4
md"## Auxiliary functions"

# ╔═╡ 45eb64c1-5af0-4987-ac1f-9d2b3dcb4c06
deltas(v::Vector)::Vector{Float64} = [0.0, [v[i] - v[i - 1] for i in 2:length(v)]...];

# ╔═╡ 1b044783-0f5f-4321-abda-35e5b7ae67c4
# R and Zsn for a given IMF, stellar metallicity and stellar yield model
function recycled_fraction(
	model::String, 
	imf::Function, 
	z::Float64,
	masses::Vector{<:Unitful.Quantity},
)::NTuple{2, Float64}

	m_low, m_high, m_ir = ustrip.(u"Msun", masses)
	sub_df = @subset(sy_data, :s_Z .== z, :model .== model)

	mass_df = @subset(sub_df, m_ir .* u"Msun" .< :s_m .< m_high .* u"Msun")
	mass = mass_df[!, "s_m"]
	
	mass_norm = [
		10 .^ range(log10(m_low), log10(m_high), step=0.1)..., m_high
	] .* u"Msun"
	
	R_int = trapz(mass, (mass .- mass_df[!, "m_rem"]) .* imf.(mass))
	Zsn_int = trapz(mass, mass .* mass_df[!, "zf_rem"] .* imf.(mass))
	norm = trapz(mass_norm, mass_norm .* imf.(mass_norm))
	
	return R_int / norm, Zsn_int / R_int
end;

# ╔═╡ 6e3dab9c-2fbe-4705-995e-753014502ede
# R and Zsn for a given IMF and stellar yield model, averaged over metallicity
function recycled_fraction(
	model::String, 
	imf::Function,
	masses::Vector{<:Unitful.Quantity},
)::NTuple{2, Float64}

	m_low, m_high, m_ir = ustrip.(u"Msun", masses)
	sub_df = @subset(sy_data, :model .== model)
	metals = unique(sub_df[!, :s_Z])
	
	R = Float64[]
	Zsn = Float64[]
	for z in metals
		mass_df = @subset(
			sub_df, 
			m_ir .* u"Msun" .< :s_m .< m_high .* u"Msun", 
			:s_Z .== z,
		)
		mass = mass_df[!, "s_m"]
		
		mass_norm = [
			10 .^ range(log10(m_low), log10(m_high), step=0.1)..., m_high
		] .* u"Msun"
		
		R_int = trapz(mass, (mass .- mass_df[!, "m_rem"]) .* imf.(mass))
		Zsn_int = trapz(mass, mass .* mass_df[!, "zf_rem"] .* imf.(mass))
		norm = trapz(mass_norm, mass_norm .* imf.(mass_norm))

		push!(R, R_int / norm)
		push!(Zsn, Zsn_int / R_int)
	end
	
	return mean(R), mean(Zsn)
end;

# ╔═╡ f5ab2c06-0d7c-4d8a-84f0-b77c97a7438d
begin
	remanent_data = DataFrame(
		imf = String[], 
		model = String[], 
		s_Z = Float64[], 
		R = Float64[], 
		Zsn = Float64[],
	)

	# [m_low, m_high_ m_ir]
	# m_low: Minimum stellar mass for the IMF population
	# m_high: Maximum stellar mass for the IMF population
	# m_ir: Lower mass limit for the instantaneous recycling hypothesis
	masses = [0.08, 100, 8] .* u"Msun"
	
	for model in values(model_names)
		sub_df = @subset(sy_data, :model .== model)
		metals = unique(sub_df[!, "s_Z"])
		for z in metals
			for imf in values(imf_funcs)
				R, Zsn = recycled_fraction(model, imf[2], z, masses)
				push!(remanent_data, [imf[1], model, z, R, Zsn])
			end
		end
	end

	sort!(remanent_data, ["imf", "model", "s_Z"])
end

# ╔═╡ 61563dc2-ee4c-4a03-9b7d-44184be7240e
begin
	MODEL_fix = model_names["Woosley1995"]
	IMF_fix = imf_funcs["Kroupa2001"][2]
	
	M_Hs = [30, 40, 50, 60, 70, 80, 90, 100]
	metal_df = @subset(sy_data, :model .== MODEL_fix)
	Zs = unique(metal_df[!, "s_Z"])

	R_mean = Float64[]
	Zsn_mean = Float64[]
	σ_R = Float64[]
	σ_Zsn = Float64[]
	
	for m_high in M_Hs
		local masses = [0.1, m_high, 8] .* u"Msun"
		remnants = Vector{Float64}[]
		for z in Zs
			R, Zsn = recycled_fraction(MODEL_fix, IMF_fix, z, masses)
			push!(remnants, [R, Zsn])
		end
		push!(R_mean, mean(getindex.(remnants, 1)))
		push!(Zsn_mean, mean(getindex.(remnants, 2)))
		push!(σ_R, std(getindex.(remnants, 1)))
		push!(σ_Zsn, std(getindex.(remnants, 2)))
	end
end

# ╔═╡ 19d23d0a-6f44-4444-8597-cd2e65e3ad1f
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"m_\mathrm{high}", 
		ylabel = L"R",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		titlesize = 25,
		xticks = M_Hs
	)

	errorbars!(ax, M_Hs, R_mean, σ_R, whiskerwidth = 10, color = :red)
	scatter!(ax, M_Hs, R_mean, markersize = 5, color = :red)

	f
end

# ╔═╡ c3717165-0646-433b-8e09-d406337a2f4c
let
	f = Figure()
	ax = Axis(
		f[1,1], 
		xlabel = L"m_\mathrm{high}", 
		ylabel = L"Z_\mathrm{SN}",
		xlabelsize = 30,
		ylabelsize = 30,
		xticklabelsize = 20,
		yticklabelsize = 20,
		titlesize = 25,
		xticks = M_Hs
	)

	errorbars!(ax, M_Hs, Zsn_mean, σ_Zsn, whiskerwidth = 10, color = :red)
	scatter!(ax, M_Hs, Zsn_mean, markersize = 5, color = :red)

	f
end

# ╔═╡ cbd6bb86-1845-4f51-bca3-59ec0e35f1af
begin  
	const NUMEQU = 5           # Number of equations
	const c₁ = ustrip(C₁)      # [Gyr Mₒ^(1/2) pc^(-3/2)]
	const c₂ = ustrip(C₂)      # [Gyr Mₒ pc^(-3)]
	const c₃ = ustrip(C₃)      # [Gyr Mₒ pc^(-3)]
	const Zeff = 1e-3 * Zsun
	const MODEL = "Woosley1995"

	# "Salpeter1955B" has a upper mass limit of 100 M⊙, which is the closest, from 
	# the Salpeter models used in PopStar, to the one referenced below
	const IMF = "Salpeter1955B"
	# From Springel et al. 2003 (https://doi.org/10.1046/j.1365-8711.2003.06206.x)
	# It uses Salpeter (1955) with slope −1.35 and upper and lower limits 
	# of 40 M⊙ and 0.1 M⊙, respectively
	const R, Zsn = recycled_fraction(
		model_names[MODEL], 
		imf_funcs["Salpeter1955"][2], 
		[0.1, 40, 8] .* u"Msun",
	) 
end

# ╔═╡ 8c7d7904-743d-44ed-bf50-b058e187b5ba
function system!(dydt, y, params, t)
	
	# Variables
	i, a, m, z, s = y
	
	# Parameters
	ρ₀, g₀, interp_ion, interp_diss = params
	
	# Auxiliary equations
	g = i + a + m
	τS = (c₁ * g₀) / sqrt(g * ρ₀)
    τR = c₂ / (i * ρ₀)
    τC = c₃ / (g * ρ₀ * (z + Zeff))
	ηdiss = interp_diss[IMF](z)
	ηion = interp_ion[IMF](z)
    recombination = i / τR
    cloud_formation = a / τC
    ψ = m / τS
	
	# ODE system
    dydt[1] = -recombination + (ηion + R) * ψ         
    dydt[2] = -cloud_formation + recombination + (ηdiss - ηion) * ψ
    dydt[3] = cloud_formation - (1 + ηdiss) * ψ
    dydt[4] = (Zsn * R - z) * ψ
	dydt[5] = (1 - R) * ψ
	
end;

# ╔═╡ bbb7263a-91e4-4a23-9e5f-416b6b7fcf6e
function integrate_model(ic, tspan, p, args = (); kwargs = (;))
	solve(ODEProblem(system!, ic, tspan, p), args...; kwargs...)
end;

# ╔═╡ eb0b2d77-b80e-43bc-b85b-86b1e38ccdd3
# System of equations adapted to compute the Jacobian
function diff_system!(dydt, y, params, t)
	
	# Variables
	i, a, m, z, s = y
	
	# Parameters
	ρ₀, g₀, ηion, ηdiss = params
	
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
	
end;

# ╔═╡ a22bc077-0e85-457b-80bc-8f9fa27b6852
# Analytical computation fo the Jacobian
begin
	@variables t y[1:NUMEQU] params[1:4]
	dydt = Vector{Num}(undef, NUMEQU)
	diff_system!(dydt, y, params, t)
	
	# Symbolic jacobian as an Array
	jac_sym_arr = Symbolics.jacobian(dydt, y)
	
	# Numerical Jacobian as an Array
	jac_num_arr = [
		build_function(jac_sym_arr[i, j], y, params, t, expression = Val{false}) for 
		i in 1:NUMEQU,
		j in 1:NUMEQU
	]
end;

# ╔═╡ 9bdbdbf7-6a1a-4bde-a12d-2181c30571bb
function jacobian!(dydt, y, parameters, t)
    for i in 1:NUMEQU
		for j in 1:NUMEQU
			dydt[(i - 1) * NUMEQU + j] = jac_num_arr[i, j](y, parameters, t)
		end
	end
end;

# ╔═╡ 1b738fec-65c4-46df-834f-55380b1db860
begin
	
	head = "int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *ode_params)
{ 
  (void)(t);

  struct ODEParameters parameters = *(struct ODEParameters *)ode_params;
  double rho0                     = parameters.rho0;
  double g0                       = parameters.g0;
  double eta_ion                  = eval_interp(parameters.interp_ion, y[3]);
  double eta_diss                 = eval_interp(parameters.interp_diss, y[3]);

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, $NUMEQU, $NUMEQU);
  gsl_matrix *m            = &dfdy_mat.matrix;\n\n"
	tail = "  dfdt[0] = 0;\n  dfdt[1] = 0;\n  dfdt[2] = 0;\n  dfdt[3] = 0;\n  dfdt[4] = 0;\n\n  return GSL_SUCCESS;\n};"

	@variables dfdt[1:(NUMEQU * NUMEQU)]
	dfdt= collect(dfdt)
	jacobian_string = ""
	
	for i in 1:NUMEQU
		for j in 1:NUMEQU
			
			index = (i - 1) * NUMEQU + j
			dfdt[index] = jac_num_arr[i, j](y, params, t)
			c_func = build_function(
				simplify(dfdt[index]), 
				y, params, t, 
				target=Symbolics.CTarget(),
			)
			
			c_func = replace(
				c_func, 
				"\n" => "",
				"RHS3" => "t",
				"RHS1" => "y",
				"RHS2[0]" => "rho0",
 				"RHS2[1]" => "g0",
  				"RHS2[2]" => "eta_ion ",
  				"RHS2[3]" => "eta_diss",
				"+ -" => "- ",
			)
			m = match(r"(?<=\{)(.*?)(?=\})", c_func)
			global jacobian_string *= replace(
				m.match,
				"  du[0] = " => "  gsl_matrix_set(m, $(i-1), $(j-1), ",
				";" => ");\n"
			)
			
		end
		global jacobian_string *= "\n"
	end

	open("./data/jacobian.txt", "w") do file
		write(file, head * jacobian_string * tail)
	end
	
end;

# ╔═╡ d51d7d61-e52f-4e42-9f52-ab31c8bf4746
function timeticks(
	t_span::Vector{<:Unitful.Quantity}, 
	n_steps::Int64,
)::Vector{Float64}

	tspan = ustrip.(u"Gyr", t_span)
	
	t_start = tspan[1] == 0.0 ? -5 : log10(tspan[1])
	t_end = log10(tspan[2])
	
	step = (t_end - t_start) / n_steps
	
	return [10^i for i in t_start:step:t_end]
	
end;

# ╔═╡ 28cfb49f-66b4-4cdf-8fcf-3f0019dff939
function phase_name_to_index(name::String)::Int64
	if name == "ionized"
		return 1
	elseif name == "atomic"
		return 2
	elseif name == "molecular"
		return 3
	elseif name == "metal"
		return 4
	elseif name == "stellar"
		return 5
	else
		@error "The phase $name does not exist!, check the spelling."
	end
end;

# ╔═╡ 8737830d-6f35-416d-bf90-cc2fbf0a4b76
# This is from the Julia source code (evalfile in base/loading.jl) but with
# the modification that it returns the module instead of the last object
function include_module(path::String)
	
	name = Symbol(basename(path))
	mod = Module(name)
	Core.eval(
		mod,
        Expr(
			:toplevel,
            :(eval(x) = $(Expr(:core, :eval))($name, x)),
            :(include(x) = $(Expr(:top, :include))($name, x)),
            :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
            :(include($path)),
		),
	)
	
	return mod
	
end;

# ╔═╡ 5d3bd48b-482b-43ed-84fb-42c2058cde74
# Compute the interpolation function for η vs. z
#
# max_age::Float64 = maximum stellar age up to which the Q values will be integrated
# ion::Bool = If the output will be the interpolation function of ηion vs. z 
# (ion = true) or ηdiss vs. z (ion = false)
function get_interp_eta(max_age::Float64, ion::Bool)::Dict
	table = Array{Float64}(undef, 6, 0)
    for imf in Q_imfs
        sol = Float64[]
        for met in Q_metals
			sub_df = @subset(imf, :Zmet .== met, $("log(age)") .< max_age)
			ages = exp10.(sub_df[!, "log(age)"]) * u"yr"
    	    Q = sub_df[!, ion ? "Q_ion" : "Q_diss"]
    	    append!(
				sol, 
				uconvert(Unitful.NoUnits, trapz(ages, Q) * (ion ? f_ion : f_diss)),
		    )
        end
        table = hcat(table, sol)
    end

	etas = DataFrame(table, col_names[2:end])

	return Dict(
		colname => LinearInterpolation(
			Q_metals, 
			ydata, 
			extrapolation_bc = Flat()
		) for (ydata, colname) in zip(eachcol(etas), col_names[2:end])
	)
end;

# ╔═╡ 97369db9-1d27-4088-b1b4-f73434ce80ea
function fractions(
	phase_name::String,
	ρ_PDF::Union{Nothing,Function}, 
	params::Params, 
	init_cond::Vector{Float64},
	t_span::Vector{<:Unitful.Quantity};
	n_steps::Int64 = 10000,
	time_ticks::Union{Nothing,Vector{Float64}} = nothing,
	log_var::Bool = true,
)::NTuple{2, Vector{Float64}}
	
	# Time values
	if time_ticks === nothing
		times = timeticks(t_span, n_steps)
	else
		times = time_ticks
	end

	# Phase to be used in the output
	phase_index = phase_name_to_index(phase_name)

	# Interpolation functions to compute ηion and ηdiss
	max_age = log10(ustrip(u"yr", t_span[2]))
	interp_eta_ion = get_interp_eta(max_age, true)
	interp_eta_diss = get_interp_eta(max_age, false)

	# Integration
	if ρ_PDF === nothing
		
		#######################################################
		# Integration using a single density value: ρ₀
		#######################################################
		
		sol = integrate_model(
			init_cond, 
			ustrip.(u"Gyr", t_span),
			[
				params.ρ₀, 
				init_cond[1] + init_cond[2] + init_cond[3], 
				interp_eta_ion, 
				interp_eta_diss,
			],
		)
		frac = sol(times)[phase_index, :]
		
	else
		
		#######################################################
		# Integration using a distribution of densities
		#######################################################

		frac = zeros(Float64, length(times))
		mass_f, points = mass_fraction(ρ_PDF, params, log_var)

		for i in 1:params.divisions	
	        sol = integrate_model(
				init_cond, 
				ustrip.(u"Gyr", t_span), 
				[
					log_var ? exp(points[i]) * params.ρ₀ : points[i] * params.ρ₀, 
					init_cond[1] + init_cond[2] + init_cond[3],
					interp_eta_ion, 
					interp_eta_diss,
				],
			)
			frac .+= sol(times)[phase_index, :] .* mass_f[i]
		end

	end	
	
	return times, frac
	
end;

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
TikzPictures = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
Trapz = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"

[compat]
CSV = "~0.10.4"
CairoMakie = "~0.8.3"
DataFrames = "~1.3.4"
DataFramesMeta = "~0.11.0"
DifferentialEquations = "~7.1.0"
HDF5 = "~0.16.9"
Interpolations = "~0.13.6"
LaTeXStrings = "~1.3.0"
LsqFit = "~0.12.1"
PlutoUI = "~0.7.39"
QuadGK = "~2.4.2"
SpecialFunctions = "~2.1.6"
Symbolics = "~4.6.0"
TikzPictures = "~3.4.2"
Trapz = "~2.0.3"
Unitful = "~1.11.0"
UnitfulAstro = "~1.1.1"

[extras]
CPUSummary = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "dd2f52bc149ff35158827471453e2e4f1a2685a6"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.26.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "7691b0012205dd7940205c54a35fd9750527d872"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.8"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d3a275e927d411e054c4192e5aca03998c233e94"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.7"

[[deps.ArrayInterfaceGPUArrays]]
deps = ["Adapt", "ArrayInterfaceCore", "GPUArrays", "LinearAlgebra"]
git-tree-sha1 = "2b6bf6e383302e0eb71e76c6fb7bd24e0e582f4d"
uuid = "6ba088a2-8465-4c0a-af30-387133b534db"
version = "0.1.0"

[[deps.ArrayInterfaceOffsetArrays]]
deps = ["ArrayInterface", "OffsetArrays", "Static"]
git-tree-sha1 = "3cbe45d8cc9cff51f302df1f87df64095423fd96"
uuid = "015c0d05-e682-4f19-8f0a-679ce4c54826"
version = "0.1.2"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "d7dc30474e73173a990eca86af76cae8790fa9f2"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.2"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c23473c60476e62579c077534b9643ec400f792b"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.6"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "960ad9a4b34380595500f60add129e178740c3a6"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.0"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "b15a6bc52594f5e4a3b825858d1089618871bf9d"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.36"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "705e7822597b432ebe152baa844b49f8026df090"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.3"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "28bbdbf0354959db89358d1d79d421ff31ef0b5e"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.3"

[[deps.BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SparseArrays"]
git-tree-sha1 = "d6a331230022493b704e1d5c11f928e2cce2d058"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.8.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "0eaf4aedad5ccc3e39481db55d72973f856dc564"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.22"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA"]
git-tree-sha1 = "5b4842a5c7e49020e25d3abe1028f8feffd636f1"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.8.3"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Chain]]
git-tree-sha1 = "339237319ef4712e6e5df7758d0bccddf5c237d9"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.10"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "c892b6671c36d64a4e3caf70ffd57f9d24f4117a"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.7"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "7297381ccb5df764549818d9a7d57e45f1057d30"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.18.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "0f4e115f6f34bbe43c19751c90a38b2f380637b9"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.3"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "87e84b2293559571802f97dd9c94cfd6be52c5e5"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.44.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DEDataArrays]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "fb2693e875ba9db2e64b684b2765e210c0d41231"
uuid = "754358af-613d-5f8d-9788-280bf1605d4c"
version = "0.2.4"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[deps.DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "OrderedCollections", "Reexport"]
git-tree-sha1 = "f1d89a07475dc4b03c08543d1c6b4b2945f33eca"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.11.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "UnPack"]
git-tree-sha1 = "078f21d61a6f43a7b3eab4620ac958183e44cee2"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.37.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterfaceCore", "ChainRulesCore", "DEDataArrays", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "2412af65ebd60b77ef0ee2cdb9f36e1e4e808d05"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.88.2"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "a442e275e9b3252dbe03db16aba435108d701f54"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.23.0"

[[deps.DiffEqJump]]
deps = ["ArrayInterfaceCore", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "08e1f6a916dafd14e731d40e6c28633a2aaae512"
uuid = "c894b116-72e5-5b58-be3c-e6d8d4ac2b12"
version = "8.5.0"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArrays", "LinearAlgebra", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "915286127c88b9306e29229cde687d46196f4060"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.11.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "28d605d9a0ac17118fe2c5e9ce0fbb76c3ceb120"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.11.0"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqJump", "DiffEqNoiseProcess", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "3f3db9365fedd5fdbecebc3cce86dfdfe5c43c50"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.1.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "0ec161f87bf4ab164ff96dfacf4be8ffff2375fd"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.62"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "d0fa82f39c2a5cdb3ee385ad52bc05c42cb4b9f0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.5"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "03b753748fd193a7f2730c02d880da27c5a24508"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExponentialUtilities]]
deps = ["ArrayInterfaceCore", "GPUArrays", "GenericSchur", "LinearAlgebra", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "2a0f80d722d7c5d94ab29c76eefbed79ded3a696"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.17.0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastBroadcast]]
deps = ["LinearAlgebra", "Polyester", "Static"]
git-tree-sha1 = "edbffa3fc9df6587927b6bc844958c61ffda8a6c"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.16"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "9267e5f50b0e12fdfd5a2455534345c4cf2c7f7a"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.14.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "4fc79c0f63ddfdcdc623a8ce36623346a7ce9ae4"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.12.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "2f18915445b248731ec5db4e4a17e451020bf21e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.30"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "b5c7fe9cea653443736d264b85466bad8c574f4a"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.9"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArrays]]
deps = ["Adapt", "LLVM", "LinearAlgebra", "Printf", "Random", "Serialization", "Statistics"]
git-tree-sha1 = "c783e8883028bf26fb05ed4022c450ef44edd875"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "8.3.2"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "4888af84657011a65afc7a564918d281612f983a"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "63561216cf6d3a5c11be429bd243f8b2b667afaf"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.7.6"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "519923f1925dbc0676f1aed54f952a194c50846f"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.2.6"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "Mmap", "Random", "Requires"]
git-tree-sha1 = "e6b1bd8339b2af5a4c2e3103f9dda65de355127e"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.16.9"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "bab67c0d1c4662d2c4be8c6007751b0b6111de5c"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.1+0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "18be5268cf415b5e27f34980ed25a7d34261aa83"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.7"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "cb7099a0109939f16a4d3b572ba8396b1f6c7c31"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.10"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "d9a03ffc2f6650bd4c831b285637929d99a4efb5"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.5"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "bcf640979ee55b652f3b01650444eb7bbe3ea837"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "336cc738f03e069ef2cac55a104eb823455dca75"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.4"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "cae5e3dfd89b209e01bcd65b3a25e74462c67ee0"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.3.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "13b16b00144816211cbf92823ded6042490eb009"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.8.1"

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "49b0c1dd5c292870577b8f58c51072bd558febb9"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.4"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "10a20c556107dc5833d3bb7c5e45c4a6e191bd28"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "4.13.0"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg", "TOML"]
git-tree-sha1 = "771bfe376249626d3ca12bcd58ba243d3f961576"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.16+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ArrayInterfaceStaticArrays", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "90d1465b8c4d22b58da9e00414f0123be3267f46"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.10.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "46a39b9c58749eefb5f2dc1178cb8fab5332b1ab"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.15"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "3da92dc6d4bce776cc850c18274a7508b8d752ef"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.7"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "f4d24f461dacac28dcd1f63ebd88a8d9d0799389"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.0"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterfaceCore", "DocStringExtensions", "GPUArrays", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SuiteSparse", "UnPack"]
git-tree-sha1 = "5812ab11c1528fef08e3d9b1fc376b1231ec8952"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.18.1"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SIMDTypes", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "af281f51681cd1b83facb8d2f1c5729c078be3e5"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.115"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "OptimBase", "Random", "StatsBase"]
git-tree-sha1 = "91aa1442e63a77f101aff01dec5a821a17f43922"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.12.1"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Serialization", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "UnicodeFun"]
git-tree-sha1 = "96e1be5153bd04212e8a9fa19b76f8eff1bb9432"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.17.3"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "cd999cfcda9ae0dd564a968087005d25359344c9"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.3.1"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test"]
git-tree-sha1 = "5c1e3d66b3a36029de4e5ac07ab8bafd5a8041e5"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.4.1"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0886d229caaa09e9f56bcf1991470bd49758a69f"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.3"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "6bb7786e4f24d44b4e29df03c69add1b63d88f01"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "393fc4d82a73c6fe0e2963dd7c882b09257be537"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "4050cd02756970414dab13b55d55ae1826b19008"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.2"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.NonlinearSolve]]
deps = ["ArrayInterfaceCore", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "8a00c7b9418270f1fa57da319d11febbe5f92101"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.20"

[[deps.Observables]]
git-tree-sha1 = "dfd8d34871bc3ad08cd16026c1828e271d554db9"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.1"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "e7fa2526bf068ad5cbfe9ba7e8a9bbd227b3211b"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.1"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "7a28efc8e34d5df89fc87343318b0a8add2c4021"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.0"

[[deps.OptimBase]]
deps = ["NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "9cb1fee807b599b5f803809e85c81b582d2009d6"
uuid = "87e2bd06-a317-5318-96d9-3ecbac512eee"
version = "2.0.2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceGPUArrays", "ArrayInterfaceStaticArrays", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "a41992efa22a08109fc8b5a4856c9c0c8a0bab90"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.14.0"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "3411935b2904d5ad3917dee58c03f0d9e6ca5355"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.11"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "e925a64b8585aa9f4e3047b8d2cdc3f0e79fd4e4"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.16"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "1155f6f937fa2b94104162f01fa400e192e4272f"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.4.2"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a121dfbba67c94a5bec9dde613c3d0cbcf3a12b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.3+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.PoissonRandom]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "44d018211a56626288b5d3f8c6497d28c26dc850"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.0"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "b5ef32913639cb417e967da27f98402938408127"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.6.11"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "7e597df97e46ffb1c8adbaddfa56908a7a20194b"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.5"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Poppler_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "e11443687ac151ac6ef6699eb75f964bed8e1faa"
uuid = "9c32591e-4766-534b-9725-b71a8799265b"
version = "0.87.0+2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "ae8a080dfc68e5615654591d8db5d53ac382353e"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.3.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "747f4261ebe38a2bc6abf0850ea8c6d9027ccd07"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "afeacaecf4ed1649555a19cb2cad3c141bbc9474"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.5.0"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArrays", "LinearAlgebra", "RecipesBase", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "23b198c159dc4fd02204683ed5d4ae105aec5053"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.29.0"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "a9a852c7ebb08e2a40e8c0ab9830a744fa283690"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.10"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "307761d71804208c0c62abdbd0ea6822aa5bbefd"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.2.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMD]]
git-tree-sha1 = "7dbc15af7ed5f751a82bf3ed37757adf76c32402"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.1"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "dd4195d308df24f33fb10dde7c22103ba88887fa"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "ac399b5b163b9140f9c310dfe9e9aaa225617ff6"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.32"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "c49318f1b9ca3d927ae576d323fa6f724d01ba53"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.1"

[[deps.SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "fe49d53d715bfbb5df78b965f5d2660781ea151c"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.38.3"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "db8481cf5d6278a121184809e9eb1628943c7704"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.13"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "f71f06cce80d21e9b93933bc04b1334bb90259ab"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.23.0"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "a9e798cae4867e3a41cae2dd9eb60c047f1212db"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.6"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "39c9f91521de844bad65049efd4f9223e7ed43f9"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.14"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "5d2c08cef80c7a3a8ba9ca023031a85c263012c5"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "cd56bf18ed715e8b09f06ef8c6b781e6cdc49911"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c82aaa13b44ea00134f8c9c89819477bd3986ecd"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.3.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "fa04638e98850332978467a085e58aababfa203a"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.8.0"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqJump", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "fea4cc29ff7d392ceb29bb64a717e6ed128bb5ff"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.49.1"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "710ec0f0e9be327b92975d7480ed6759ade69ec3"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.3.9"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "9abba8f8fb8458e9adf07c8a2377a070674a24f1"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.8"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "Reexport", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "6549d3b1b5cf86446949c62616675588159ea2fb"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.9.4"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "04777432d74ec5bc91ca047c9e0e0fd7f81acdb6"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.1+0"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "bfa211c9543f8c062143f2a48e5bcbb226fd790b"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.7"

[[deps.Symbolics]]
deps = ["ArrayInterfaceCore", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "90872c5192844bc6b25fcba733b7b0020a509258"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.6.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Tectonic]]
deps = ["Pkg"]
git-tree-sha1 = "0b3881685ddb3ab066159b2ce294dc54fcf3b9ee"
uuid = "9ac5f52a-99c6-489f-af81-462ef484790f"
version = "0.8.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "d223de97c948636a4f34d1f84d92fd7602dc555b"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.10"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "f90022b44b7bf97952756a6b6737d1a0024a3233"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.5"

[[deps.TikzPictures]]
deps = ["LaTeXStrings", "Poppler_jll", "Requires", "Tectonic"]
git-tree-sha1 = "4e75374d207fefb21105074100034236fceed7cb"
uuid = "37f6aa50-8035-52d0-81c2-5a1d08754b2d"
version = "3.4.2"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "7638550aaea1c9a1e86817a231ef0faa9aca79bd"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.19"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c76399a3bbe6f5a88faa33c8f8a65aa631d95013"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.73"

[[deps.Trapz]]
git-tree-sha1 = "79eb0ed763084a3e7de81fe1838379ac6a23b6a0"
uuid = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
version = "2.0.3"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "b8d08f55b02625770c09615d96927b3a8396925e"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.11"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b649200e887a487468b71821e2644382699f1b0f"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.11.0"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "d6cfdb6ddeb388af1aea38d2b9905fa014d92d98"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.6.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "c4e1c470a94063b911fd1b1a204cd2bb34a8cd15"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.1.1"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "c81cd227862cf6892d4751c7842aaa856f10269a"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.35"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╠═982068e0-59bb-11ec-27f5-51126c2ba1df
# ╟─08df960b-fd82-43ba-a9dc-bf5e83af587e
# ╟─3b1726c6-60c2-45be-932d-efa8d2ef23e0
# ╟─6173713c-92b0-43ec-8713-1cbf442aa1ce
# ╟─a842b24e-8d26-41ab-9de3-91632aede893
# ╟─64787011-b5b8-42be-b6e4-37ebc5138b3e
# ╟─76fe97bd-36c8-40d2-9b5a-0ea5059bd7c7
# ╟─14c7f574-0623-4254-b8f7-97984d32351c
# ╟─bc9ab101-7cc3-4baa-b83d-ce546f6b576d
# ╟─047bbd39-9cf9-4bd7-b38e-16aa505b0b08
# ╟─2fe0dc4c-da44-4fc8-bef8-1fa615a0fe4a
# ╟─eaf272c7-4162-4a9a-92e3-9835c6158394
# ╟─ac553b12-4857-4cc1-8ea2-fe9e8863b429
# ╟─dc6fd12b-c821-4e20-a896-25c8aab9df94
# ╟─1d27ec35-65ca-4c94-9e8d-54d1c11e759f
# ╟─68732d91-805a-4663-9166-f8483213a8d2
# ╟─27281e53-e519-4ad0-af5d-59fb0e208534
# ╟─6503fb74-c34f-40db-afb4-7affd4ceef88
# ╟─f6251e55-f88b-4f53-8449-e30b0bf9ae44
# ╟─7099a821-a0f0-4931-b6cf-88581e9cff9e
# ╟─4a7eb24b-0874-49a3-9b08-4ffb6a7f0ce7
# ╟─f2a6676f-457a-476a-9ce7-c336aa9bf47f
# ╟─040e1a8c-97ab-4751-a556-ed936fe58c35
# ╟─3767c7f9-a0bc-467a-a20a-5e5a266111c7
# ╟─127e1dfa-62d8-4721-adc8-cb24c6e9cdcc
# ╟─005957d6-6f27-4abc-a872-45cf6a032b9f
# ╟─0fcd2ad5-440c-4128-be21-1f8a354074fe
# ╟─f8b02d00-ff30-480e-b5eb-e150e4678c95
# ╟─44c88ad8-a8c3-45e3-9a56-be3ce5bf66fa
# ╟─448e1dee-4628-4c14-9d6f-dc165b2e826e
# ╟─c7409abf-dc22-429e-ad4d-e2cbd465d454
# ╟─aa5e9990-db35-4a91-912e-f839daf6c686
# ╟─342c1ad8-338e-44e2-adec-7638fe1767a2
# ╟─7788b98a-5bec-4b6d-82d9-2c272e2255a7
# ╟─7af2d4a6-a304-404d-8cbe-eeddb80beba6
# ╟─aeb72f0e-2252-486a-b79b-9d8cc6e5f962
# ╟─e83337bd-8c2d-4a9a-bd8b-7f8201cf67ad
# ╟─b8beaaa6-b018-4ca8-b19d-a918dc761707
# ╟─43eb3af9-86c5-49e9-af0e-3270e3df493e
# ╟─05d22c2d-f76f-4931-8e21-6d31e9ab178e
# ╟─ea0bada1-4359-4c2b-9dc8-d91b6ebd5686
# ╟─00e1f4f0-6a6a-48dc-81d5-99f356aa3410
# ╟─0f212d88-6b3a-4f96-8f05-d5b0f9943fe3
# ╟─3c26c7f7-0d0c-4fe3-a071-ab76c7904659
# ╟─fa38bb0f-c807-48b2-8bf2-33e457b53435
# ╟─e879e854-a9c3-4762-9963-dfa8960f2dc5
# ╟─66f233f2-76d8-45be-9aa2-831ff7269d96
# ╟─d7069cdb-9a99-4194-93ea-03b5b82ddc89
# ╟─4525cc1d-5338-400e-840a-47c069014ec6
# ╟─0e350f3a-534e-4b7b-b473-9f3ebefb399e
# ╟─fd2f6a01-9125-49fb-8615-b2e5deccd107
# ╟─930fe232-f4ff-4331-9752-4e4ec0e66009
# ╟─386334a4-60bf-47de-bf72-516f042b1407
# ╟─9b0f3856-bee3-485e-9049-1a082c86e571
# ╟─d5ba04de-e2e8-44af-9b60-6a47b782248e
# ╟─396aa4fe-5d65-4852-9cb2-03c654201f6e
# ╟─399a13c4-235b-40e7-8a2f-5affd889c014
# ╟─f6c88e38-5f80-4ab7-afef-4f3249af8723
# ╟─ff8a6018-dd86-4b09-a122-a72e0dfa7013
# ╟─2c58b32c-e731-467b-b051-7063b3d3e341
# ╟─2b7d684e-3db3-4966-9993-0446a7db7edb
# ╟─0a465f30-6904-4195-830e-21cfb00fb63a
# ╟─d7871e31-d841-4951-aec5-d6f2915d0ceb
# ╟─9b738a7d-1a85-460b-9585-6ee4f32b41ae
# ╟─afdb1375-ac94-4dcc-8fcb-c732e6029b83
# ╟─1e2f1eff-3295-4eae-ac83-569e9f8211fd
# ╟─67dfb1b3-4818-445f-bfe2-16d442506567
# ╟─3da12a3d-bc26-4115-b769-78841499e434
# ╟─64fda1f0-8fcf-4898-b3dd-883f151f69b3
# ╟─76136e33-14f7-4c5e-84e9-970372c6c01a
# ╟─1fce7934-0ce4-4059-ba8e-691f1c505f4d
# ╟─e9eca322-1466-4d18-af62-28d3edd42cfc
# ╟─68502224-d186-464b-8d35-1ce9ad0a9994
# ╟─f017b433-05dd-48bc-b0b4-f70a39100b2d
# ╟─533b3cd0-c1f6-4ecd-b196-4ed35bf77135
# ╟─97e7f37b-9494-4ae9-a076-77b90a974a81
# ╟─be85ba3b-5439-4cf3-bb14-d24d61a283c3
# ╟─d1e89b59-bc6f-46fd-a7cd-126fad530916
# ╟─f5ab2c06-0d7c-4d8a-84f0-b77c97a7438d
# ╟─49d39360-3609-407c-bfee-c46e7485727a
# ╟─61563dc2-ee4c-4a03-9b7d-44184be7240e
# ╟─19d23d0a-6f44-4444-8597-cd2e65e3ad1f
# ╟─c3717165-0646-433b-8e09-d406337a2f4c
# ╟─9666bdc8-cbc0-4757-9bd8-a76477c252eb
# ╟─ca9a233b-d3ca-4a76-a3d8-f29884ac9484
# ╠═cbd6bb86-1845-4f51-bca3-59ec0e35f1af
# ╟─e2e4ae4f-dcdc-4999-88f2-853378be859a
# ╟─2be98f2a-57a2-4f53-ad53-b1c0e9e9aafa
# ╠═8c7d7904-743d-44ed-bf50-b058e187b5ba
# ╟─01c985a1-4d3d-4799-8e9e-5e25390a9fbb
# ╠═eb0b2d77-b80e-43bc-b85b-86b1e38ccdd3
# ╠═a22bc077-0e85-457b-80bc-8f9fa27b6852
# ╠═9bdbdbf7-6a1a-4bde-a12d-2181c30571bb
# ╠═1b738fec-65c4-46df-834f-55380b1db860
# ╟─ceec9b81-6a43-4bb9-bb74-b309ef4c3037
# ╟─4607856c-7472-4131-a2ee-29f7150f5cb4
# ╠═bbb7263a-91e4-4a23-9e5f-416b6b7fcf6e
# ╠═7a2987ef-d37e-4c7a-aaa8-8186694bea88
# ╠═97369db9-1d27-4088-b1b4-f73434ce80ea
# ╟─4cfe1c80-c67e-4dd3-825b-d893800d68c0
# ╟─d7ba9e0c-5cfa-4176-adff-12cb8e20679b
# ╠═82e78dc9-b89e-48d9-9f70-6f3238dfd196
# ╟─34c053bf-0b4f-45c4-bb79-7e5e89a26060
# ╠═768f8ffb-a08b-4498-97e4-1a3a866e69c7
# ╟─8864b4d3-6a9f-4e7b-8cd6-ed32a0116f4a
# ╠═aad5e227-67a5-49a0-a79d-24160d3ebe06
# ╠═d6502bbe-ea58-471d-a748-d9c8f37ce61e
# ╟─b3a260b6-eb31-43a0-9fd6-60a507984319
# ╠═5ba3a0c1-6107-45a1-9b1d-5c323b9a7145
# ╟─946d007d-abd7-4cd3-9789-e77b1ad6ebf4
# ╠═45eb64c1-5af0-4987-ac1f-9d2b3dcb4c06
# ╠═1b044783-0f5f-4321-abda-35e5b7ae67c4
# ╠═6e3dab9c-2fbe-4705-995e-753014502ede
# ╠═d51d7d61-e52f-4e42-9f52-ab31c8bf4746
# ╠═28cfb49f-66b4-4cdf-8fcf-3f0019dff939
# ╠═8737830d-6f35-416d-bf90-cc2fbf0a4b76
# ╠═5d3bd48b-482b-43ed-84fb-42c2058cde74
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

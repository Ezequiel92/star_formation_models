/***************************************************************************************
* Numerical solution for the full model (variable τS without infall and outflow).
****************************************************************************************/

const char name[] = "full_model_variable_tau_no_IO";

vector_type Y0{5}; // Initial condition vector.

/* Constants. */

const Double K = 1 / 400.0;		  // [Mₒ^(-2/3) pc^(4/3) Gyr^(-1)]
const Double Twarm = 10000.0;	  // [K]
const Double Tcold = 100.0;		  // [K]
const Double T1 = 50000.0;		  // [K]
const Double Sigma_ion = 8e-4;	  // [Mₒ pc^(-2)]
const Double Sigma_diss = 1.5e-4; // [Mₒ pc^(-2)]
const Double C2 = 0.074;		  // [Mₒ^2 pc^(-4) Gyr]
const Double C4 = 798e-3;		  // [Mₒ^2 pc^(-4) Gyr]
const Double R = 0.17;
const Double alpha = 7.0;
const Double eta_i_lim = 1000.0;
const Double eta_d_lim = 100.0;
const Double Zsun = 0.006;
const Double Zsn = 0.06;
const Double Zeff = Zsun * 1e-3;

/* Star formation rate function */

Double SFR(Double i, Double a, Double m, Double s)
{
	Double star_elem = m + alpha * a;
	Double g = i + a + m;
	Double tot = g + s;
	Double tau_S = pow(star_elem, 1 / 3.0) / (pow(g, 1 / 3.0) * pow(tot, 2 / 3.0) * K);

	return (star_elem / tau_S); // [Mₒ pc^(-2) Gyr^(-1)]
}

/* Witness function. */

struct file_observer
{
	FILE *fp;

	file_observer(FILE *file) : fp{file} {}

	template <class State>
	void operator()(const State &y, Double t) const
	{	// t	i(t)	a(t)	m(t)	z(t)	SFR(t)	s(t)	g(t) = i + a + m
        fprintf(fp, "%.15g \t %.15g \t %.15g \t %.15g \t %.15g \t %.15g \t %.15g \t %.15g\n", 
				t, y[0], y[1], y[2], y[4], SFR(y[0], y[1], y[2], y[3]), y[3], y[0] + y[1] + y[2]);
	}
};

/* Equations. */

struct rhs
{
	template <class State>
	void operator()(const State &y, State &dydt, Double t)
	{
		/*
		*	ionized gas:       i(t) -> y[0]
		*	atomic gas:        a(t) -> y[1]
		*	molecular gas:     m(t) -> y[2]
		*	star mass:         s(t) -> y[3]
		*	metals:            z(t) -> y[4]

		*	Units -> Each equation has Mₒ pc^(-2) Gyr^(-1) [Solar_mass * parsec^-2 * years^-9)] as units in the LHS and RHS.
		*/

		(void)(t);

		Double star_elem = y[2] + alpha * y[1];
		Double g = y[0] + y[1] + y[2];
		Double tot = g + y[3];
		Double tau_S = pow(star_elem, 1 / 3.0) / (pow(g, 1 / 3.0) * pow(tot, 2 / 3.0) * K);
		Double psi = star_elem / tau_S;
		Double eta_ion = eta_i_lim * (1 - exp(-y[1] / Sigma_ion));
		Double eta_diss = eta_d_lim * (1 - exp(-y[2] / Sigma_diss));
		Double Z = y[4] / g;
		Double tau_R = C2 * (1 + T1 * psi / (Twarm * g)) / (g * tot);
		Double tau_C = C4 * (1 + T1 * psi / (Tcold * g)) * Zsun / (g * tot * (Z + Zeff));
		Double recombination = y[0] / tau_R;
		Double cloud_formation = y[1] / tau_C;

		dydt[0] = -recombination + (eta_ion + R) * psi;
		dydt[1] = -cloud_formation + recombination + (eta_diss - eta_ion - alpha * y[1] / star_elem) * psi;
		dydt[2] = cloud_formation - (eta_diss + y[2] / star_elem) * psi;
		dydt[3] = (1 - R) * psi;
		dydt[4] = (Zsn * R - Z) * psi;
	}
};

/* Jacobian. */

struct jacobi
{
	template <class State, class Matrix>
	void operator()(const State &y, Matrix &J, const Double &t, State &dfdt)
	{
		(void)(t);
		J(0, 0) = (0.00166666666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.000833333333333333 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (13.5135135135135 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (13.5135135135135 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (13.5135135135135 * y[0] * (y[1] + y[0] + y[2] + y[3])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (13.5135135135135 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(0, 1) = (0.00166666666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.000833333333333333 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.0116666666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) + (3125. * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / exp(1250. * y[1]) + (13.5135135135135 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.0583333333333333 * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666)))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (13.5135135135135 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (13.5135135135135 * y[0] * (y[1] + y[0] + y[2] + y[3])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(0, 2) = (0.00166666666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.000833333333333333 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.00166666666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) + (13.5135135135135 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.00833333333333333 * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666)))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (13.5135135135135 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (13.5135135135135 * y[0] * (y[1] + y[0] + y[2] + y[3])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(0, 3) = (0.00166666666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.112612612612613 * y[0] * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (13.5135135135135 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(0, 4) = 0;

		J(1, 0) = (0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.000833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (13.5135135135135 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (13.5135135135135 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (13.5135135135135 * y[0] * (y[1] + y[0] + y[2] + y[3])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (13.5135135135135 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * y[4]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667)) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(1, 1) = (0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + 0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1.25e6 / exp(1250. * y[1]) + (49. * y[1]) / pow(7. * y[1] + y[2], 2) - 7. / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) + (0.000833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.0116666666666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) - (13.5135135135135 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.0583333333333333 * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666)))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (13.5135135135135 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (13.5135135135135 * y[0] * (y[1] + y[0] + y[2] + y[3])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * y[4]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (5.83333333333333 * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666))) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (208.855472013367 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(1, 2) = (0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + 0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (666666.666666667 / exp(6666.66666666666667 * y[2]) + (7. * y[1]) / pow(7. * y[1] + y[2], 2)) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) + (0.000833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.00166666666666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) - (13.5135135135135 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.00833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.00833333333333333 * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666)))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (13.5135135135135 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (13.5135135135135 * y[0] * (y[1] + y[0] + y[2] + y[3])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * y[4]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.833333333333333 * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666))) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(1, 3) = (0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) - (0.112612612612613 * y[0] * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (13.5135135135135 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (174.046226677806 * y[1] * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(1, 4) = (-208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3])) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));

		J(2, 0) = (-0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) - (0.000833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (100. * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * y[4]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667)) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(2, 1) = (-0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.0175 * y[2] * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(7. * y[1] + y[2], 1.3333333333333333) - (0.000833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (100. * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (0.0116666666666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * y[4]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (5.83333333333333 * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666))) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472013367 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(2, 2) = (-0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) - 0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (666666.666666667 / exp(6666.66666666666667 * y[2]) - (1. * y[2]) / pow(7. * y[1] + y[2], 2) + 1 / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) - (0.000833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (100. * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (0.00166666666666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * y[4]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) - (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[3]) * ((0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333)) - (0.833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.833333333333333 * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666))) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(2, 3) = (-0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) - (174.046226677806 * y[1] * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (208.855472013367 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[4] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
		J(2, 4) = (208.855472013367 * y[1] * (y[1] + y[0] + y[2] + y[3])) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));

		J(3, 0) = (0.00138333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.000691666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666);
		J(3, 1) = (0.00138333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.000691666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.00968333333333333 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333);
		J(3, 2) = (0.00138333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.000691666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.00138333333333333 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333);
		J(3, 3) = (0.00138333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333);
		J(3, 4) = 0;

		J(4, 0) = (0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) * y[4]) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0102 - (1. * y[4]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.000833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) * (0.0102 - (1. * y[4]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2], 0.6666666666666666);
		J(4, 1) = (0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) * y[4]) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0102 - (1. * y[4]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.000833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) * (0.0102 - (1. * y[4]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.0116666666666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) * (0.0102 - (1. * y[4]) / (y[1] + y[0] + y[2]))) / pow(7. * y[1] + y[2], 0.3333333333333333);
		J(4, 2) = (0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) * y[4]) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0102 - (1. * y[4]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333) + (0.000833333333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) * (0.0102 - (1. * y[4]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.00166666666666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666) * (0.0102 - (1. * y[4]) / (y[1] + y[0] + y[2]))) / pow(7. * y[1] + y[2], 0.3333333333333333);
		J(4, 3) = (0.00166666666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0102 - (1. * y[4]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[3], 0.3333333333333333);
		J(4, 4) = (-0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[3], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666);

		dfdt[0] = 0;
		dfdt[1] = 0;
		dfdt[2] = 0;
		dfdt[3] = 0;
		dfdt[4] = 0;
	}
};
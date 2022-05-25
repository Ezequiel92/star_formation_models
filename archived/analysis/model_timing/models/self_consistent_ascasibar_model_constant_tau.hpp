/***************************************************************************************
* Numerical solution for the Self consistent Ascasibar et al. model (constant τS).
****************************************************************************************/

const char name[] = "Self consistent Ascasibar et al. model";

/* Constants. */

const Double tau_S = 2.6;  // [Gyr]
const Double C1 = 3.03e-6; // [Gyr]
const Double C2 = 8.74e-5; // [Gyr]
const Double R = 0.18;
const Double eta_ion = 955.29;
const Double eta_diss = 380.93;
const Double Zsun = 0.0134;
const Double Zsn = 0.09;
const Double Zeff = Zsun * 1e-3;

/* Initial condition for test run and time measurement. */

const Double g0 = 1.0;
const Double i0 = g0 * 0.6;	 // [Mₒ pc^(-2)]
const Double a0 = g0 * 0.2;	 // [Mₒ pc^(-2)]
const Double m0 = g0 * 0.2;	 // [Mₒ pc^(-2)]
const Double s0 = g0 * 0.0;	 // [Mₒ pc^(-2)]
const Double z0 = g0 * 1e-4; // [Mₒ pc^(-2)]
vector_type Y0{5};

/* Parameters. */

const Double Tstart = 0.0; // [Gyr]
const Double Tend = 1.0;   // [Gyr]

/* Star formation rate function */

Double SFR(Double m)
{
	return m / tau_S; // [Mₒ pc^(-3) Gyr^(-1)]
}

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

		*	Units -> Each equation has Mₒ pc^(-3) Gyr^(-1) [Solar_mass * parsec^-3 * years^-9)] as units in the LHS and RHS.
		*/

		(void)(t);

		Double g = y[0] + y[1] + y[2];
		Double psi = y[2] / tau_S;
		Double Z = y[4] / g;
		Double tau_R = C1 / y[0];
		Double tau_C = C2 / ((y[1] + y[2]) * (Z + Zeff));
		Double recombination = y[0] / tau_R;
		Double cloud_formation = y[1] / tau_C;

		dydt[0] = -recombination + (eta_ion + R) * psi;
		dydt[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
		dydt[2] = cloud_formation - (eta_diss + 1) * psi;
		dydt[3] = (1 - R) * psi;
		dydt[4] = (Zsn * R - Z) * psi;
	}
};

/* Constant Jacobian. */

struct jacobi_const
{
	template <class State, class Matrix>
	void operator()(const State &y, Matrix &J, const Double &t, State &dfdt)
	{
		(void)(y);
		(void)(t);

		J(0, 0) = -6.6006600660066e6;
		J(0, 1) = 0;
		J(0, 2) = 367.488461538462;
		J(0, 3) = 0;
		J(0, 4) = 0;

		J(1, 0) = 6.60177013465649e6;
		J(1, 1) = -1445.97254004577;
		J(1, 2) = -382.372223200141;
		J(1, 3) = 0;
		J(1, 4) = -1.11006864988558e7;

		J(2, 0) = -1110.06864988558;
		J(2, 1) = 1445.97254004577;
		J(2, 2) = 14.5683770462947;
		J(2, 3) = 0;
		J(2, 4) = 1.11006864988558e7;

		J(3, 0) = 0;
		J(3, 1) = 0;
		J(3, 2) = 0.315384615384615;
		J(3, 3) = 0;
		J(3, 4) = 0;

		J(4, 0) = 3.84615384615385e-7;
		J(4, 1) = 3.84615384615385e-7;
		J(4, 2) = 0.00619269230769231;
		J(4, 3) = 0;
		J(4, 4) = -0.00384615384615385;

		dfdt[0] = 0;
		dfdt[1] = 0;
		dfdt[2] = 0;
		dfdt[3] = 0;
		dfdt[4] = 0;
	}
};

/* Variable Jacobian. */

struct jacobi_var
{
	template <class State, class Matrix>
	void operator()(const State &y, Matrix &J, const Double &t, State &dfdt)
	{
		(void)(t);

		J(0, 0) = -660066.00660066 * y[0];
		J(0, 1) = 0;
		J(0, 2) = 367.488461538462;
		J(0, 3) = 0;
		J(0, 4) = 0;

		J(1, 0) = 660066.00660066 * y[0] + (11441.647597254 * y[1] * (y[1] + y[2]) * y[4]) / pow(y[1] + y[0] + y[2], 2);
		J(1, 1) = (11441.647597254 * y[1] * (y[1] + y[2]) * y[4]) / pow(y[1] + y[0] + y[2], 2) - 11441.647597254 * y[1] * (0.0000134 + y[4] / (y[1] + y[0] + y[2])) - 11441.647597254 * (y[1] + y[2]) * (0.0000134 + y[4] / (y[1] + y[0] + y[2]));
		J(1, 2) = -220.907692307692 + (11441.647597254 * y[1] * (y[1] + y[2]) * y[4]) / pow(y[1] + y[0] + y[2], 2) - 11441.647597254 * y[1] * (0.0000134 + y[4] / (y[1] + y[0] + y[2]));
		J(1, 3) = 0;
		J(1, 4) = (-11441.647597254 * y[1] * (y[1] + y[2])) / (y[1] + y[0] + y[2]);

		J(2, 0) = (-11441.647597254 * y[1] * (y[1] + y[2]) * y[4]) / pow(y[1] + y[0] + y[2], 2);
		J(2, 1) = (-11441.647597254 * y[1] * (y[1] + y[2]) * y[4]) / pow(y[1] + y[0] + y[2], 2) + 11441.647597254 * y[1] * (0.0000134 + y[4] / (y[1] + y[0] + y[2])) + 11441.647597254 * (y[1] + y[2]) * (0.0000134 + y[4] / (y[1] + y[0] + y[2]));
		J(2, 2) = -146.896153846154 - (11441.647597254 * y[1] * (y[1] + y[2]) * y[4]) / pow(y[1] + y[0] + y[2], 2) + 11441.647597254 * y[1] * (0.0000134 + y[4] / (y[1] + y[0] + y[2]));
		J(2, 3) = 0;
		J(2, 4) = (11441.647597254 * y[1] * (y[1] + y[2])) / (y[1] + y[0] + y[2]);

		J(3, 0) = 0;
		J(3, 1) = 0;
		J(3, 2) = 0.315384615384615;
		J(3, 3) = 0;
		J(3, 4) = 0;

		J(4, 0) = (0.384615384615385 * y[2] * y[4]) / pow(y[1] + y[0] + y[2], 2);
		J(4, 1) = (0.384615384615385 * y[2] * y[4]) / pow(y[1] + y[0] + y[2], 2);
		J(4, 2) = (0.384615384615385 * y[2] * y[4]) / pow(y[1] + y[0] + y[2], 2) + 0.384615384615385 * (0.0162 - (1. * y[4]) / (y[1] + y[0] + y[2]));
		J(4, 3) = 0;
		J(4, 4) = (-0.384615384615385 * y[2]) / (y[1] + y[0] + y[2]);

		dfdt[0] = 0;
		dfdt[1] = 0;
		dfdt[2] = 0;
		dfdt[3] = 0;
		dfdt[4] = 0;
	}
};
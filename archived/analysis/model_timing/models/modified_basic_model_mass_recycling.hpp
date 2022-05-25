/***************************************************************************************
* Numerical solution for the modified basic model with mass recycling.
* (modified to have the same constants as the Ascasibar model)
****************************************************************************************/

const char name[] = "Modified basic model with mass recycling";

/* Constants. */

const Double tau_S = 2.6;       // [Gyr]
const Double n = 1e-3;          // [1/cm^3]
const Double Rh = 1.9;          // [cm^3 Gyr^-1]
const Double SigmaNuB = 8158.0; // [cm^3 Gyr^-1]
const Double Zsun = 0.0134;
const Double Zeff = 1e-3 * Zsun;
const Double Zsn = 0.09;
const Double eta_ion = 955.29;
const Double eta_diss = 380.93;
const Double R = 0.18;

/* Initial condition for test run and time measurement. */

const Double i0 = 0.6;
const Double a0 = 0.2;
const Double m0 = 0.2;
const Double s0 = 0.0;
const Double z0 = 1e-4;
vector_type Y0{5};

/* Parameters. */

const Double Z = z0;
const Double Tstart = 0.0; // [Gyr]
const Double Tend = 10.0;   // [Gyr]

/* Star formation rate function (fractional -> SFR / g) */

Double SFRf(Double mf)
{
    return (mf / tau_S); // [Gyr^(-1)]
}

/* Equations. */

struct rhs
{
    template <class State>
    void operator()(const State &y, State &dydt, Double t)
    {
        /*
		*	ionized gas fraction:       i(t)/g -> y[0]
		*	atomic gas fraction:        a(t)/g -> y[1]
		*	molecular gas fraction:     m(t)/g -> y[2]
		*	star fraction:              s(t)/g -> y[3]
		*	metal fraction:             z(t)/g -> y[4]

		*	Units -> Each equation has Gyr^(-1) [years^-9] as units in the LHS and RHS.
		*/

        (void)(t);

		Double ne = n;
		Double nh = n;
        Double recombination = y[0] * ne * SigmaNuB;
        Double cloud_formation = y[1] * 2 * nh * Rh * (Z + Zeff) / Zsun;
        Double psi = y[2] / tau_S;

        dydt[0] = -recombination + eta_ion * psi + R * psi;
        dydt[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
        dydt[2] = cloud_formation - (1 + eta_diss) * psi;
		dydt[3] = psi - R * psi;
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

		J(0,0) = -8.158;
		J(0,1) = 0;
		J(0,2) = 367.41923076923075;
		J(0,3) = 0;
		J(0,4) = 0;

		J(1,0) = 8.158;
		J(1,1) = -0.0000321582089552239;
		J(1,2) = -220.9076923076923;
		J(1,3) = 0;
		J(1,4) = 0;

		J(2,0) = 0;
		J(2,1) = 0.0000321582089552239;
		J(2,2) = -146.89615384615382;
		J(2,3) = 0;
		J(2,4) = 0;

		J(3,0) = 0;
		J(3,1) = 0;
		J(3,2) = 0.3846153846153846;
		J(3,3) = 0;
		J(3,4) = 0;

		J(4,0) = 0;
		J(4,1) = 0;
		J(4,2) = 0.006192307692307691;
		J(4,3) = 0;
		J(4,4) = 0;

		dfdt[0] = 0;
		dfdt[1] = 0;
		dfdt[2] = 0;
		dfdt[3] = 0;
		dfdt[4] = 0;
    }
};


/*  Only for compatibility with the solver API. */

struct jacobi_var
{
    template <class State, class Matrix>
    void operator()(const State &y, Matrix &J, const Double &t, State &dfdt)
    {
        (void)(y);
        (void)(t);
        (void)(J);
        (void)(dfdt);
    }
};
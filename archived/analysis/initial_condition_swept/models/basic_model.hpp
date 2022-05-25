/***************************************************************************************
* Numerical solution for the basic model.
****************************************************************************************/

const char name[] = "basic_model";

vector_type Y0{4}; // Initial condition vector.

/* Constants. */

const Double tau_S = 3.0;       // [Gyr]
const Double Rh = 1.9;          // [cm^3 Gyr^-1]
const Double SigmaNuB = 8158.0; // [cm^3 Gyr^-1]
const Double Zsun = 0.02;
const Double Zeff = 1e-3 * Zsun;
const Double Zsn = 0.2;
const Double eta_ion = 955.0;
const Double eta_diss = 381.0;
const Double R = 0.17;

/* Star formation rate function (fractional -> SFR / g) */

Double SFRf(Double mf)
{
    return (mf / tau_S); // [Gyr^(-1)]
}

/* Witness function. */

struct file_observer
{
    FILE *fp;

    file_observer(FILE *file) : fp{file} {}

    template <class State>
    void operator()(const State &y, Double t) const
    {	// t	if(t)	af(t)	mf(t)	zf(t)	SFRf(t)
        fprintf(fp, "%.15g \t %.15g \t %.15g \t %.15g \t %.15g \t %.15g \n", 
				t, y[0], y[1], y[2], y[3], SFRf(y[2]));
    }
};

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
		*	metal fraction:             z(t)/g -> y[3]

		*	Units -> Each equation has Gyr^(-1) [years^-9] as units in the LHS and RHS.
		*/

        (void)(t);

        Double recombination = y[0] * n * SigmaNuB;
        Double cloud_formation = y[1] * 2 * n * Rh * (Z + Zeff) / Zsun;
        Double psi = y[2] / tau_S;

        dydt[0] = -recombination + eta_ion * psi;
        dydt[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
        dydt[2] = cloud_formation - (1 + eta_diss) * psi;
        dydt[3] = (Zsn * R - Z) * psi;
    }
};

/* Jacobian. */

struct jacobi
{
    template <class State, class Matrix>
    void operator()(const State &y, Matrix &J, const Double &t, State &dfdt)
    {
        (void)(y);
        (void)(t);

        J(0, 0) = -8.158;
        J(0, 1) = 0;
        J(0, 2) = 318.333333333333;
        J(0, 3) = 0;

        J(1, 0) = 8.158;
        J(1, 1) = -0.0000228;
        J(1, 2) = -191.333333333333;
        J(1, 3) = 0;

        J(2, 0) = 0;
        J(2, 1) = 0.0000228;
        J(2, 2) = -127.333333333333;
        J(2, 3) = 0;

        J(3, 0) = 0;
        J(3, 1) = 0;
        J(3, 2) = 0.0113;
        J(3, 3) = 0;

        dfdt[0] = 0;
        dfdt[1] = 0;
        dfdt[2] = 0;
        dfdt[3] = 0;
    }
};

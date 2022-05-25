/***************************************************************************************
* Numerical solution for the simplified model.
****************************************************************************************/

const char name[] = "Simplified model";

/* Constants. */

const Double K = 1 / 400.0;       // [Mₒ^(-2/3) pc^(4/3) Gyr^(-1)]
const Double Twarm = 10000.0;     // [K]
const Double Tcold = 100.0;       // [K]
const Double T1 = 50000.0;        // [K]
const Double C2 = 0.074;          // [Mₒ^2 pc^(-4) Gyr]
const Double C4 = 798e-3;         // [Mₒ^2 pc^(-4) Gyr]
const Double Sigma_ion = 8e-4;    // [Mₒ pc^(-2)]
const Double Sigma_diss = 1.5e-4; // [Mₒ pc^(-2)]
const Double R = 0.17;
const Double alpha = 7.0;
const Double eta_i_lim = 955.29;
const Double eta_d_lim = 380.93;
const Double Zsun = 0.0134;
const Double Zsn = 0.2;
const Double Zeff = Zsun * 1e-3;

/* Initial condition for test run and time measurement. */

const Double i0 = 0.6;  // [Mₒ pc^(-2)]
const Double a0 = 0.2;  // [Mₒ pc^(-2)]
const Double m0 = 0.2;  // [Mₒ pc^(-2)]
const Double s0 = 0.0;  // [Mₒ pc^(-2)]
const Double z0 = 1e-4; // [Mₒ pc^(-2)]
vector_type Y0{5};

/* Parameters. */

const Double g = 1.0;      // [Mₒ pc^(-2)]
const Double Tstart = 0.0; // [Gyr]
const Double Tend = 10.0;   // [Gyr]

/* File initialization function. */

FILE *file_ini(const char *name)
{
    FILE *file;
    file = fopen(name, "a");
    fprintf(file, "t \t if(t) \t af(t) \t mf(t) \t sf(t) \t zf(t) \t SFRf(t) \n");
    return file;
}

/* Star formation rate function (fractional -> SFR / g) */

Double SFRf(Double af, Double mf)
{
    Double star_elem = mf + alpha * af;

    return (K * pow(star_elem * g, 2 / 3.0)); // [Gyr^(-1)]
}

/* Witness function. */

struct file_observer
{
    FILE *fp;

    file_observer(FILE *file) : fp{file} {}

    template <class State>
    void operator()(const State &y, Double t) const
    {
        fprintf(fp, "%.15g \t %.15g \t %.15g \t %.15g \t %.15g \t %.15g \t %.15g \n", t, y[0], y[1], y[2], y[3], y[4], SFRf(y[1], y[2]));
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
		*	star fraction:              s(t)/g -> y[3]
		*	metal fraction:             z(t)/g -> y[4]

		*	Units -> Each equation has Gyr^(-1) [years^-9] as units in the LHS and RHS.
		*/

        (void)(t);

        Double star_elem = y[2] + alpha * y[1];
        Double psi_f = K * pow(star_elem * g, 2 / 3.0);
        Double Z = y[4];
        Double eta_ion = eta_i_lim * (1 - exp(-g * y[1] / Sigma_ion));
        Double eta_diss = eta_d_lim * (1 - exp(-g * y[2] / Sigma_diss));
        Double tau_R = C2 * (1 + T1 * psi_f / Twarm) / (g * g);
        Double tau_C = C4 * (1 + T1 * psi_f / Tcold) * Zsun / (g * g * (Z + Zeff));
        Double recombination = y[0] / tau_R;
        Double cloud_formation = y[1] / tau_C;

        dydt[0] = -recombination + (eta_ion + R) * psi_f;
        dydt[1] = -cloud_formation + recombination + (eta_diss - eta_ion - alpha * y[1] / star_elem) * psi_f;
        dydt[2] = cloud_formation - (eta_diss + y[2] / star_elem) * psi_f;
        dydt[3] = (1 - R) * psi_f;
        dydt[4] = (Zsn * R - Z) * psi_f;
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

        J(0, 0) = -13.2863206266003;
        J(0, 1) = 9.92148263978142;
        J(0, 2) = 1.41735466282592;
        J(0, 3) = 0;
        J(0, 4) = 0;

        J(1, 0) = 13.2863206266003;
        J(1, 1) = -6.13313477783284;
        J(1, 2) = -0.873465602465963;
        J(1, 3) = 0;
        J(1, 4) = -6.90171258848667;

        J(2, 0) = 0;
        J(2, 1) = -3.79662699549041;
        J(2, 2) = -0.545071793723071;
        J(2, 3) = 0;
        J(2, 4) = 6.90171258848667;

        J(3, 0) = 0;
        J(3, 1) = 0.00827913354182634;
        J(3, 2) = 0.00118273336311805;
        J(3, 3) = 0;
        J(3, 4) = 0;

        J(4, 0) = 0;
        J(4, 1) = 0.000338147743455317;
        J(4, 2) = 0.0000483068204936167;
        J(4, 3) = 0;
        J(4, 4) = -0.00341995189335339;

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
        J(0, 0) = -13.5135135135135 / (1. + 0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666));
        J(0, 1) = (0.0116666666666667 * (0.17 + 955.29 * (1. - 1. / exp(1250. * y[1])))) / pow(7. * y[1] + y[2], 0.3333333333333333) + (2985.28125 * pow(7. * y[1] + y[2], 0.6666666666666666)) / exp(1250. * y[1]) + (0.788288288288288 * y[0]) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(1. + 0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666), 2));
        J(0, 2) = (0.00166666666666667 * (0.17 + 955.29 * (1. - 1. / exp(1250. * y[1])))) / pow(7. * y[1] + y[2], 0.3333333333333333) + (0.112612612612613 * y[0]) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(1. + 0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666), 2));
        J(0, 3) = 0;
        J(0, 4) = 0;

        J(1, 0) = 13.5135135135135 / (1. + 0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666));
        J(1, 1) = 0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * (-1.1941125e6 / exp(1250. * y[1]) + (49. * y[1]) / pow(7. * y[1] + y[2], 2) - 7. / (7. * y[1] + y[2])) + (0.0116666666666667 * (-955.29 * (1. - 1. / exp(1250. * y[1])) + 380.93 * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(7. * y[1] + y[2], 0.3333333333333333) - (0.788288288288288 * y[0]) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(1. + 0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666), 2)) + (545.51802391551 * y[1] * (0.0000134 + y[4])) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(1. + 1.25 * pow(7. * y[1] + y[2], 0.6666666666666666), 2)) - (93.5173755283732 * (0.0000134 + y[4])) / (1. + 1.25 * pow(7. * y[1] + y[2], 0.6666666666666666));
        J(1, 2) = 0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * (2.53953333333333e6 / exp(6666.66666666666667 * y[2]) + (7. * y[1]) / pow(7. * y[1] + y[2], 2)) + (0.00166666666666667 * (-955.29 * (1. - 1. / exp(1250. * y[1])) + 380.93 * (1. - 1. / exp(6666.66666666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(7. * y[1] + y[2], 0.3333333333333333) - (0.112612612612613 * y[0]) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(1. + 0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666), 2)) + (77.9311462736443 * y[1] * (0.0000134 + y[4])) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(1. + 1.25 * pow(7. * y[1] + y[2], 0.6666666666666666), 2));
        J(1, 3) = 0;
        J(1, 4) = (-93.5173755283732 * y[1]) / (1. + 1.25 * pow(7. * y[1] + y[2], 0.6666666666666666));

        J(2, 0) = 0;
        J(2, 1) = (0.0175 * y[2]) / pow(7. * y[1] + y[2], 1.3333333333333333) - (0.0116666666666667 * (380.93 * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(7. * y[1] + y[2], 0.3333333333333333) - (545.51802391551 * y[1] * (0.0000134 + y[4])) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(1. + 1.25 * pow(7. * y[1] + y[2], 0.6666666666666666), 2)) + (93.5173755283732 * (0.0000134 + y[4])) / (1. + 1.25 * pow(7. * y[1] + y[2], 0.6666666666666666));
        J(2, 2) = -0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * (2.53953333333333e6 / exp(6666.66666666666667 * y[2]) - (1. * y[2]) / pow(7. * y[1] + y[2], 2) + 1 / (7. * y[1] + y[2])) - (0.00166666666666667 * (380.93 * (1. - 1. / exp(6666.66666666666667 * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(7. * y[1] + y[2], 0.3333333333333333) - (77.9311462736443 * y[1] * (0.0000134 + y[4])) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(1. + 1.25 * pow(7. * y[1] + y[2], 0.6666666666666666), 2));
        J(2, 3) = 0;
        J(2, 4) = (93.5173755283732 * y[1]) / (1. + 1.25 * pow(7. * y[1] + y[2], 0.6666666666666666));

        J(3, 0) = 0;
        J(3, 1) = 0.00968333333333333 / pow(7. * y[1] + y[2], 0.3333333333333333);
        J(3, 2) = 0.00138333333333333 / pow(7. * y[1] + y[2], 0.3333333333333333);
        J(3, 3) = 0;
        J(3, 4) = 0;

        J(4, 0) = 0;
        J(4, 1) = (0.0116666666666667 * (0.034 - 1. * y[4])) / pow(7. * y[1] + y[2], 0.3333333333333333);
        J(4, 2) = (0.00166666666666667 * (0.034 - 1. * y[4])) / pow(7. * y[1] + y[2], 0.3333333333333333);
        J(4, 3) = 0;
        J(4, 4) = -0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666);

        dfdt[0] = 0;
        dfdt[1] = 0;
        dfdt[2] = 0;
        dfdt[3] = 0;
        dfdt[4] = 0;
    }
};
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/*******************************************************************************************
* Model 014
*
* Alternative names: 'Simplified model'
*******************************************************************************************/

static const char NAME[] = "Model 014"; // Model name
static const int SIZE = 5;              // Size of the system of ODEs

/* File initialization function */

FILE *file_ini(const char *name)
{
    FILE *file;
    file = fopen(name, "a");
    fprintf(file, "t \t i(t) \t a(t) \t m(t) \t z(t) \t s(t) \n");
    return file;
}

/* Constants. */

const double K = 1 / 400.0;       // [Mₒ^(-2/3) pc^(4/3) Gyr^(-1)]
const double Twarm = 10000.0;     // [K]
const double Tcold = 100.0;       // [K]
const double T1 = 50000.0;        // [K]
const double Sigma_ion = 8e-4;    // [Mₒ pc^(-2)]
const double Sigma_diss = 1.5e-4; // [Mₒ pc^(-2)]
const double C2 = 0.074;          // [Mₒ^2 pc^(-4) Gyr]
const double C4 = 798e-3;         // [Mₒ^2 pc^(-4) Gyr]
const double alpha = 7.0;
const double eta_i_lim = 955.29;
const double eta_d_lim = 380.93;
const double Zsun = 0.0134;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.2;
const double R = 0.17;

/* System of equations */

/*
	Ionized gas fraction:       i(t) / g -> y[0]
	Atomic gas fraction:        a(t) / g -> y[1]
	Molecular gas fraction:     m(t) / g -> y[2]
	Metal fraction:             z(t) / g -> y[3]
	Stellar fraction:           s(t) / g -> y[4]		

    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
*/
static int ode_system(double t, const double y[], double f[], void *params)
{
    (void)(t);

    /* Parameters */
    double g = *(double *)params;

    /* Auxiliary equations */
    double star_elem = y[2] + alpha * y[1];
    double psi = K * pow(star_elem * g, 2 / 3.0);
    double eta_ion = eta_i_lim * (1 - exp(-g * y[1] / Sigma_ion));
    double eta_diss = eta_d_lim * (1 - exp(-g * y[2] / Sigma_diss));
    double Z = y[3];
    double tau_R = C2 * (1 + T1 * psi / Twarm) / (g * g);
    double tau_C = C4 * (1 + T1 * psi / Tcold) * Zsun / (g * g * (Z + Zeff));
    double recombination = y[0] / tau_R;
    double cloud_formation = y[1] / tau_C;

    /* ODE system */
    f[0] = -recombination + (eta_ion + R) * psi;
    f[1] = -cloud_formation + recombination + (eta_diss - eta_ion - alpha * y[1] / star_elem) * psi;
    f[2] = cloud_formation - (eta_diss + y[2] / star_elem) * psi;
    f[3] = (Zsn * R - Z) * psi;
    f[4] = (1 - R) * psi;

    return GSL_SUCCESS;
};

/* Jacobian */

int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    (void)(t);

    double g = *(double *)params;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 5, 5);
    gsl_matrix *m = &dfdy_mat.matrix;

    gsl_matrix_set(m, 0, 0, (-13.51351351 * pow(g, 2)) / (1. + 0.0125 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666)));
    gsl_matrix_set(m, 0, 1, (0.01166666667 * (0.17 + 955.29 * (1. - 1. / exp(1250. * g * y[1]))) * g) / pow(g * (7. * y[1] + y[2]), 0.3333333333333333) + (2985.28125 * g * pow(g * (7. * y[1] + y[2]), 0.6666666666666666)) / exp(1250. * g * y[1]) + (0.7882882883 * pow(g, 3) * y[0]) / (pow(g * (7. * y[1] + y[2]), 0.3333333333333333) * pow(1. + 0.0125 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666), 2)));
    gsl_matrix_set(m, 0, 2, (0.001666666667 * (0.17 + 955.29 * (1. - 1. / exp(1250. * g * y[1]))) * g) / pow(g * (7. * y[1] + y[2]), 0.3333333333333333) + (0.1126126126 * pow(g, 3) * y[0]) / (pow(g * (7. * y[1] + y[2]), 0.3333333333333333) * pow(1. + 0.0125 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666), 2)));
    gsl_matrix_set(m, 0, 3, 0);
    gsl_matrix_set(m, 0, 4, 0);

    gsl_matrix_set(m, 1, 0, (13.51351351 * pow(g, 2)) / (1. + 0.0125 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666)));
    gsl_matrix_set(m, 1, 1, 0.0025 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666) * ((-1.1941125e6 * g) / exp(1250. * g * y[1]) + (49. * y[1]) / pow(7. * y[1] + y[2], 2) - 7. / (7. * y[1] + y[2])) + (0.01166666667 * g * (-955.29 * (1. - 1. / exp(1250. * g * y[1])) + 380.93 * (1. - 1. / exp(6666.666666667 * g * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(g * (7. * y[1] + y[2]), 0.3333333333333333) - (0.7882882883 * pow(g, 3) * y[0]) / (pow(g * (7. * y[1] + y[2]), 0.3333333333333333) * pow(1. + 0.0125 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666), 2)) + (545.5180239 * pow(g, 3) * y[1] * (0.0000134 + y[3])) / (pow(g * (7. * y[1] + y[2]), 0.3333333333333333) * pow(1. + 1.25 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666), 2)) - (93.51737553 * pow(g, 2) * (0.0000134 + y[3])) / (1. + 1.25 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666)));
    gsl_matrix_set(m, 1, 2, 0.0025 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666) * ((2.539533333e6 * g) / exp(6666.666666667 * g * y[2]) + (7. * y[1]) / pow(7. * y[1] + y[2], 2)) + (0.001666666667 * g * (-955.29 * (1. - 1. / exp(1250. * g * y[1])) + 380.93 * (1. - 1. / exp(6666.666666667 * g * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(g * (7. * y[1] + y[2]), 0.3333333333333333) - (0.1126126126 * pow(g, 3) * y[0]) / (pow(g * (7. * y[1] + y[2]), 0.3333333333333333) * pow(1. + 0.0125 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666), 2)) + (77.93114627 * pow(g, 3) * y[1] * (0.0000134 + y[3])) / (pow(g * (7. * y[1] + y[2]), 0.3333333333333333) * pow(1. + 1.25 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666), 2)));
    gsl_matrix_set(m, 1, 3, (-93.51737553 * pow(g, 2) * y[1]) / (1. + 1.25 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666)));
    gsl_matrix_set(m, 1, 4, 0);

    gsl_matrix_set(m, 2, 0, 0);
    gsl_matrix_set(m, 2, 1, (0.0175 * y[2] * pow(g * (7. * y[1] + y[2]), 0.6666666666666666)) / pow(7. * y[1] + y[2], 2) - (0.01166666667 * g * (380.93 * (1. - 1. / exp(6666.666666667 * g * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(g * (7. * y[1] + y[2]), 0.3333333333333333) - (545.5180239 * pow(g, 3) * y[1] * (0.0000134 + y[3])) / (pow(g * (7. * y[1] + y[2]), 0.3333333333333333) * pow(1. + 1.25 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666), 2)) + (93.51737553 * pow(g, 2) * (0.0000134 + y[3])) / (1. + 1.25 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666)));
    gsl_matrix_set(m, 2, 2, -0.0025 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666) * ((2.539533333e6 * g) / exp(6666.666666667 * g * y[2]) - (1. * y[2]) / pow(7. * y[1] + y[2], 2) + 1 / (7. * y[1] + y[2])) - (0.001666666667 * g * (380.93 * (1. - 1. / exp(6666.666666667 * g * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(g * (7. * y[1] + y[2]), 0.3333333333333333) - (77.93114627 * pow(g, 3) * y[1] * (0.0000134 + y[3])) / (pow(g * (7. * y[1] + y[2]), 0.3333333333333333) * pow(1. + 1.25 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666), 2)));
    gsl_matrix_set(m, 2, 3, (93.51737553 * pow(g, 2) * y[1]) / (1. + 1.25 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666)));
    gsl_matrix_set(m, 2, 4, 0);

    gsl_matrix_set(m, 3, 0, 0);
    gsl_matrix_set(m, 3, 1, (0.01166666667 * g * (0.034 - 1. * y[3])) / pow(g * (7. * y[1] + y[2]), 0.3333333333333333));
    gsl_matrix_set(m, 3, 2, (0.001666666667 * g * (0.034 - 1. * y[3])) / pow(g * (7. * y[1] + y[2]), 0.3333333333333333));
    gsl_matrix_set(m, 3, 3, -0.0025 * pow(g * (7. * y[1] + y[2]), 0.6666666666666666));
    gsl_matrix_set(m, 3, 4, 0);

    gsl_matrix_set(m, 4, 0, 0);
    gsl_matrix_set(m, 4, 1, (0.009683333333 * g) / pow(g * (7. * y[1] + y[2]), 0.3333333333333333));
    gsl_matrix_set(m, 4, 2, (0.001383333333 * g) / pow(g * (7. * y[1] + y[2]), 0.3333333333333333));
    gsl_matrix_set(m, 4, 3, 0);
    gsl_matrix_set(m, 4, 4, 0);

    dfdt[0] = 0;
    dfdt[1] = 0;
    dfdt[2] = 0;
    dfdt[3] = 0;
    dfdt[4] = 0;

    return GSL_SUCCESS;
};

/*******************************************************************************************
* Solving the system
*******************************************************************************************/

/* Integration span */

double t_start = 0.0; // [Gyr]
double t_end = 1.0;   // [Gyr]

int main()
{
    /* Initial time step */
    double dt = (t_end - t_start) * 1e-4;

    /* Absolute tolerance */
    static const double abs_tol = 1.0e-6;

    /* Relative tolerance */
    static const double rel_tol = 1.0e-8;

    /* Initial conditions */
    const double i0 = 0.6;
    const double a0 = 0.2;
    const double m0 = 0.2;
    const double z0 = 1e-4;
    const double s0 = 0.0;
    double Y0[] = {i0, a0, m0, z0, s0};

    /* Parameters */
    double params[] = {1.0}; // {g [Mₒ pc^(-2)]}

    /* Integration */
    gsl_odeiv2_system sys = {ode_system, jacobian, SIZE, params};
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
        &sys,
        gsl_odeiv2_step_msbdf,
        dt,
        abs_tol,
        rel_tol);

    /* Save data in file */

    FILE *file;
    file = file_ini("plots/c_msbdf.dat");
    int status;

    for (int i = 0; i <= (int)(1e4); ++i)
    {
        status = gsl_odeiv2_driver_apply(driver, &t_start, i * dt, Y0);

        if (status != GSL_SUCCESS)
        {
            printf("error, return value = %d\n", status);
            break;
        }

        fprintf(
            file,
            "%.15g \t %.15g \t %.15g \t %.15g \t %.15g \t %.15g \n",
            t_start, Y0[0], Y0[1], Y0[2], Y0[3], Y0[4]);
    }

    gsl_odeiv2_driver_free(driver);
    fclose(file);

    printf("Stellar fraction after %gGyr (msbdf): %.3g", t_end, Y0[4]);
}

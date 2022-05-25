#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/*******************************************************************************************
* Model 013
*
* Alternative names: 'Full model (Kennicutt-Schmidt law without infall and outflow)'
*******************************************************************************************/

static const char NAME[] = "Model 013"; // Model name
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

const double A = 2.5e-4; // [Mₒ pc^(-2) Gyr^(-1)]
const double N = 1.4;
const double Twarm = 10000.0;     // [K]
const double Tcold = 100.0;       // [K]
const double T1 = 50000.0;        // [K]
const double Sigma_ion = 8e-4;    // [Mₒ pc^(-2)]
const double Sigma_diss = 1.5e-4; // [Mₒ pc^(-2)]
const double C2 = 0.074;          // [Mₒ^2 pc^(-4) Gyr]
const double C4 = 798e-3;         // [Mₒ^2 pc^(-4) Gyr]
const double alpha = 0.3;
const double eta_i_lim = 955.29;
const double eta_d_lim = 380.93;
const double Zsun = 0.0134;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.2;
const double R = 0.17;

/* System of equations */

/*
	Ionized gas density:       i(t) -> y[0]
	Atomic gas density:        a(t) -> y[1]
	Molecular gas density:     m(t) -> y[2]
	Metal density:             z(t) -> y[3]
	Stellar density:           s(t) -> y[4]		

    Each equation has Mₒ pc^(-2) Gyr^(-1) [solar_mass * parsec^(-2) * years^(-9)]
	as its units on the LHS and RHS
*/
static int ode_system(double t, const double y[], double f[], void *params)
{
    /* Auxiliary equations */
    double g = y[0] + y[1] + y[2];
    double tot = g + y[4];
    double star_elem = y[2] + alpha * y[1];
    double psi = A * pow((y[1] + y[2]), N);
    double eta_ion = eta_i_lim * (1 - exp(-y[1] / Sigma_ion));
    double eta_diss = eta_d_lim * (1 - exp(-y[2] / Sigma_diss));
    double Z = y[3] / g;
    double tau_R = C2 * (1 + T1 * psi / (Twarm * g)) / (g * tot);
    double tau_C = C4 * (1 + T1 * psi / (Tcold * g)) * Zsun / (g * tot * (Z + Zeff));
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
    (void)(params);

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 5, 5);
    gsl_matrix *m = &dfdy_mat.matrix;

    gsl_matrix_set(m, 0, 0, (-13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) - (0.01689189189 * y[0] * pow(y[1] + y[2], 1.4) * (y[1] + y[0] + y[2] + y[4])) / ((y[1] + y[0] + y[2]) * pow(1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2)) - (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) - (13.51351351 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 0, 1, 0.00035 * (0.17 + 955.29 * (1. - 1. / exp(1250. * y[1]))) * pow(y[1] + y[2], 0.4) + (298.528125 * pow(y[1] + y[2], 1.4)) / exp(1250. * y[1]) - (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (13.51351351 * y[0] * (y[1] + y[0] + y[2]) * ((-0.00125 * pow(y[1] + y[2], 1.4)) / pow(y[1] + y[0] + y[2], 2) + (0.00175 * pow(y[1] + y[2], 0.4)) / (y[1] + y[0] + y[2])) * (y[1] + y[0] + y[2] + y[4])) / pow(1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2) - (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 0, 2, 0.00035 * (0.17 + 955.29 * (1. - 1. / exp(1250. * y[1]))) * pow(y[1] + y[2], 0.4) - (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (13.51351351 * y[0] * (y[1] + y[0] + y[2]) * ((-0.00125 * pow(y[1] + y[2], 1.4)) / pow(y[1] + y[0] + y[2], 2) + (0.00175 * pow(y[1] + y[2], 0.4)) / (y[1] + y[0] + y[2])) * (y[1] + y[0] + y[2] + y[4])) / pow(1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2) - (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 0, 3, 0);
    gsl_matrix_set(m, 0, 4, (-13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));

    gsl_matrix_set(m, 1, 0, (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (0.01689189189 * y[0] * pow(y[1] + y[2], 1.4) * (y[1] + y[0] + y[2] + y[4])) / ((y[1] + y[0] + y[2]) * pow(1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2)) + (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (13.51351351 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]))) - (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) - (11.68967194 * y[1] * pow(y[1] + y[2], 1.4) * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / ((y[1] + y[0] + y[2]) * pow(1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2)) - (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 1, 1, 0.00025 * pow(y[1] + y[2], 1.4) * (-1.1941125e6 / exp(1250. * y[1]) + (0.09 * y[1]) / pow(0.3 * y[1] + y[2], 2) - 0.3 / (0.3 * y[1] + y[2])) + 0.00035 * pow(y[1] + y[2], 0.4) * (-955.29 * (1. - 1. / exp(1250. * y[1])) + 380.93 * (1. - 1. / exp(6666.666666667 * y[2])) - (0.3 * y[1]) / (0.3 * y[1] + y[2])) + (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) - (13.51351351 * y[0] * (y[1] + y[0] + y[2]) * ((-0.00125 * pow(y[1] + y[2], 1.4)) / pow(y[1] + y[0] + y[2], 2) + (0.00175 * pow(y[1] + y[2], 0.4)) / (y[1] + y[0] + y[2])) * (y[1] + y[0] + y[2] + y[4])) / pow(1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2) + (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]))) - (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * ((-0.125 * pow(y[1] + y[2], 1.4)) / pow(y[1] + y[0] + y[2], 2) + (0.175 * pow(y[1] + y[2], 0.4)) / (y[1] + y[0] + y[2])) * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2) - (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) - (93.51737553 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 1, 2, 0.00025 * pow(y[1] + y[2], 1.4) * (2.539533333e6 / exp(6666.666666667 * y[2]) + (0.3 * y[1]) / pow(0.3 * y[1] + y[2], 2)) + 0.00035 * pow(y[1] + y[2], 0.4) * (-955.29 * (1. - 1. / exp(1250. * y[1])) + 380.93 * (1. - 1. / exp(6666.666666667 * y[2])) - (0.3 * y[1]) / (0.3 * y[1] + y[2])) + (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) - (13.51351351 * y[0] * (y[1] + y[0] + y[2]) * ((-0.00125 * pow(y[1] + y[2], 1.4)) / pow(y[1] + y[0] + y[2], 2) + (0.00175 * pow(y[1] + y[2], 0.4)) / (y[1] + y[0] + y[2])) * (y[1] + y[0] + y[2] + y[4])) / pow(1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2) + (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]))) - (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * ((-0.125 * pow(y[1] + y[2], 1.4)) / pow(y[1] + y[0] + y[2], 2) + (0.175 * pow(y[1] + y[2], 0.4)) / (y[1] + y[0] + y[2])) * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2) - (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 1, 3, (-93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 1, 4, (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.00125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) - (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));

    gsl_matrix_set(m, 2, 0, (-93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]))) + (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (11.68967194 * y[1] * pow(y[1] + y[2], 1.4) * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / ((y[1] + y[0] + y[2]) * pow(1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2)) + (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 2, 1, (0.000075 * y[2] * pow(y[1] + y[2], 1.4)) / pow(0.3 * y[1] + y[2], 2) - 0.00035 * pow(y[1] + y[2], 0.4) * (380.93 * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (0.3 * y[1] + y[2])) - (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]))) + (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) - (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * ((-0.125 * pow(y[1] + y[2], 1.4)) / pow(y[1] + y[0] + y[2], 2) + (0.175 * pow(y[1] + y[2], 0.4)) / (y[1] + y[0] + y[2])) * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2) + (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) + (93.51737553 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 2, 2, -0.00025 * pow(y[1] + y[2], 1.4) * (2.539533333e6 / exp(6666.666666667 * y[2]) - (1. * y[2]) / pow(0.3 * y[1] + y[2], 2) + 1 / (0.3 * y[1] + y[2])) - 0.00035 * pow(y[1] + y[2], 0.4) * (380.93 * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (0.3 * y[1] + y[2])) - (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]))) + (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])) - (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * ((-0.125 * pow(y[1] + y[2], 1.4)) / pow(y[1] + y[0] + y[2], 2) + (0.175 * pow(y[1] + y[2], 0.4)) / (y[1] + y[0] + y[2])) * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]), 2) + (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 2, 3, (93.51737553 * y[1] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 2, 4, (93.51737553 * y[1] * (y[1] + y[0] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (0.125 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2])));

    gsl_matrix_set(m, 3, 0, (0.00025 * pow(y[1] + y[2], 1.4) * y[3]) / pow(y[1] + y[0] + y[2], 2));
    gsl_matrix_set(m, 3, 1, (0.00025 * pow(y[1] + y[2], 1.4) * y[3]) / pow(y[1] + y[0] + y[2], 2) + 0.00035 * pow(y[1] + y[2], 0.4) * (0.034 - (1. * y[3]) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 3, 2, (0.00025 * pow(y[1] + y[2], 1.4) * y[3]) / pow(y[1] + y[0] + y[2], 2) + 0.00035 * pow(y[1] + y[2], 0.4) * (0.034 - (1. * y[3]) / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 3, 3, (-0.00025 * pow(y[1] + y[2], 1.4)) / (y[1] + y[0] + y[2]));
    gsl_matrix_set(m, 3, 4, 0);

    gsl_matrix_set(m, 4, 0, 0);
    gsl_matrix_set(m, 4, 1, 0.0002905 * pow(y[1] + y[2], 0.4));
    gsl_matrix_set(m, 4, 2, 0.0002905 * pow(y[1] + y[2], 0.4));
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
    const double g0 = 1e-3 * 14.8285; // [Mₒ pc^(-2)]
    const double i0 = g0 * 0.6;       // [Mₒ pc^(-2)]
    const double a0 = g0 * 0.2;       // [Mₒ pc^(-2)]
    const double m0 = g0 * 0.2;       // [Mₒ pc^(-2)]
    const double z0 = g0 * 1e-4;      // [Mₒ pc^(-2)]
    const double s0 = g0 * 0.0;       // [Mₒ pc^(-2)]
    double Y0[] = {i0, a0, m0, z0, s0};

    /* Parameters */

    /* 
	*	There are no parameters in this model.
	*/

    void *params = NULL;

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

    printf("Stellar density after %gGyr (msbdf): %.3g Mₒpc^(-2)", t_end, Y0[4]);
}

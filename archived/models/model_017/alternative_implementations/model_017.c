#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/*******************************************************************************************
* Model 017
*
* Alternative names: 'New full model', 'Model D'
*******************************************************************************************/

static const char NAME[] = "Model 017"; // Model name
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

const double K = 1 / 400.0; // [Mₒ^(-2/3) pc^(4/3) Gyr^(-1)]
const double C1 = 3.03e-6;  // [Gyr]
const double C2 = 8.74e-5;  // [Gyr]
const double eta_ion = 955.29;
const double eta_diss = 380.93;
const double Zsun = 0.0134;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.09;
const double R = 0.18;

/* System of equations */

/*
	Ionized gas density:       i(t) -> y[0]
	Atomic gas density:        a(t) -> y[1]
	Molecular gas density:     m(t) -> y[2]
	Metal density:             z(t) -> y[3]
	Stellar density:           s(t) -> y[4]		

    Each equation has Mₒ pc^(-3) Gyr^(-1) [solar_mass * parsec^(-3) * years^(-9)]
	as its units on the LHS and RHS
*/
static int ode_system(double t, const double y[], double f[], void *params)
{
    /* Auxiliary equations */
    double g = y[0] + y[1] + y[2];
    double tot = g + y[4];
    double tau_S = pow(y[2], 1 / 3.0) / (K * pow(g, 1 / 3.0) * pow(tot, 2 / 3.0));
    double psi = y[2] / tau_S;
    double Z = y[3] / g;
    double tau_R = C1 / y[0];
    double tau_C = C2 / ((y[1] + y[2]) * (Z + Zeff));
    double recombination = y[0] / tau_R;
    double cloud_formation = y[1] / tau_C;

    /* ODE system */
    f[0] = -recombination + (eta_ion + R) * psi;
    f[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
    f[2] = cloud_formation - (eta_diss + 1) * psi;
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

    gsl_matrix_set(m, 0, 0, -660066.0066 * y[0] + (1.59245 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.796225 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
    gsl_matrix_set(m, 0, 1, (1.59245 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.796225 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
    gsl_matrix_set(m, 0, 2, (1.59245 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.796225 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (1.59245 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[2], 0.3333333333333333));
    gsl_matrix_set(m, 0, 3, 0);
    gsl_matrix_set(m, 0, 4, (1.59245 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333));

    gsl_matrix_set(m, 1, 0, 660066.0066 * y[0] - (0.9572666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) - (0.4786333333 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2));
    gsl_matrix_set(m, 1, 1, (-0.9572666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) - (0.4786333333 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2) - 11441.6476 * y[1] * (0.0000134 + y[3] / (y[1] + y[0] + y[2])) - 11441.6476 * (y[1] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 1, 2, (-0.9572666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) - (0.4786333333 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (0.9572666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[2], 0.3333333333333333) + (11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2) - 11441.6476 * y[1] * (0.0000134 + y[3] / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 1, 3, (-11441.6476 * y[1] * (y[1] + y[2])) / (y[1] + y[0] + y[2]));
    gsl_matrix_set(m, 1, 4, (-0.9572666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333));

    gsl_matrix_set(m, 2, 0, (-0.63655 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) - (0.318275 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2));
    gsl_matrix_set(m, 2, 1, (-0.63655 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) - (0.318275 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2) + 11441.6476 * y[1] * (0.0000134 + y[3] / (y[1] + y[0] + y[2])) + 11441.6476 * (y[1] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 2, 2, (-0.63655 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) - (0.318275 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (0.63655 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[2], 0.3333333333333333) - (11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2) + 11441.6476 * y[1] * (0.0000134 + y[3] / (y[1] + y[0] + y[2])));
    gsl_matrix_set(m, 2, 3, (11441.6476 * y[1] * (y[1] + y[2])) / (y[1] + y[0] + y[2]));
    gsl_matrix_set(m, 2, 4, (-0.63655 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333));

    gsl_matrix_set(m, 3, 0, (0.0025 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * y[3]) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.001666666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0162 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0008333333333 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (0.0162 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
    gsl_matrix_set(m, 3, 1, (0.0025 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * y[3]) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.001666666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0162 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0008333333333 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (0.0162 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
    gsl_matrix_set(m, 3, 2, (0.0025 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * y[3]) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.001666666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0162 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0008333333333 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (0.0162 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.001666666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (0.0162 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[2], 0.3333333333333333));
    gsl_matrix_set(m, 3, 3, (-0.0025 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
    gsl_matrix_set(m, 3, 4, (0.001666666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0162 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333));

    gsl_matrix_set(m, 4, 0, (0.001366666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0006833333333 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
    gsl_matrix_set(m, 4, 1, (0.001366666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0006833333333 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
    gsl_matrix_set(m, 4, 2, (0.001366666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0006833333333 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.001366666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[2], 0.3333333333333333));
    gsl_matrix_set(m, 4, 3, 0);
    gsl_matrix_set(m, 4, 4, (0.001366666667 * pow(y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333));

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
    const double g0 = 1e-3 * 14.8285; // [Mₒ pc^(-3)]
    const double i0 = g0 * 0.6;       // [Mₒ pc^(-3)]
    const double a0 = g0 * 0.2;       // [Mₒ pc^(-3)]
    const double m0 = g0 * 0.2;       // [Mₒ pc^(-3)]
    const double z0 = g0 * 1e-4;      // [Mₒ pc^(-3)]
    const double s0 = g0 * 0.0;       // [Mₒ pc^(-3)]
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

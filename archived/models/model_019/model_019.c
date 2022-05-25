#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/*******************************************************************************************
* Model 019
*******************************************************************************************/

static const char NAME[] = "Model 019"; // Model name
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

const double C1 = 0.809;    // [Gyr Mₒ^(1/2) pc^(-3/2)]
const double C2 = 3.024e-6; // [Gyr Mₒ pc^(-3)]
const double C3 = 8.745e-5; // [Gyr Mₒ pc^(-3)]
const double Zsun = 0.0134;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.09;
const double eta_ion = 955.29;
const double eta_diss = 380.93;
const double R = 0.18;

/* System of equations */

/*
    Ionized gas fraction:       i(t) / ρ -> y[0]
    Atomic gas fraction:        a(t) / ρ -> y[1]
    Molecular gas fraction:     m(t) / ρ -> y[2]
    Metal fraction:             z(t) / ρ -> y[3]
    Stellar fraction:           s(t) / ρ -> y[4]	

    where ρ = i(t) + a(t) + m(t) + s(t)
     
    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
*/
static int ode_system(double t, const double y[], double f[], void *params)
{
    (void)(t);

    double rho0 = *(double *)params;
    double g0 = *((double *)params + 1);

    /* Auxiliary equations */
    double g = y[0] + y[1] + y[2];
    double tau_S = (C1 * g0) / sqrt(g * rho0);
    double tau_R = C2 / (y[0] * rho0);
    double tau_C = C3 / (g * rho0 * (y[3] + Zeff));
    double recombination = y[0] / tau_R;
    double cloud_formation = y[1] / tau_C;
    double psi = y[2] / tau_S;

    /* ODE system */
    f[0] = -recombination + (eta_ion + R) * psi;
    f[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
    f[2] = cloud_formation - (1 + eta_diss) * psi;
    f[3] = (Zsn * R - y[3]) * psi;
    f[4] = (1 - R) * psi;

    return GSL_SUCCESS;
};

/* Jacobian */

int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    (void)(t);

    double rho0 = *(double *)params;
    double g0 = *((double *)params + 1);

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 5, 5);
    gsl_matrix *m = &dfdy_mat.matrix;

    gsl_matrix_set(m, 0, 0, -6.61e5 * rho0 * y[0] + (590. * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))));
    gsl_matrix_set(m, 0, 1, (590. * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))));
    gsl_matrix_set(m, 0, 2, (590. * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) + (1180. * sqrt(rho0 * (y[1] + y[0] + y[2]))) / g0);
    gsl_matrix_set(m, 0, 3, 0);
    gsl_matrix_set(m, 0, 4, 0);

    gsl_matrix_set(m, 1, 0, 6.61e5 * rho0 * y[0] - (354.7 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) - 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]));
    gsl_matrix_set(m, 1, 1, (-354.7 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) - 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]) - 1.143e4 * rho0 * (y[1] + y[0] + y[2]) * (0.0000134 + y[3]));
    gsl_matrix_set(m, 1, 2, (-354.7 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) - (709.3 * sqrt(rho0 * (y[1] + y[0] + y[2]))) / g0 - 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]));
    gsl_matrix_set(m, 1, 3, -1.143e4 * rho0 * y[1] * (y[1] + y[0] + y[2]));
    gsl_matrix_set(m, 1, 4, 0);

    gsl_matrix_set(m, 2, 0, (-235.8 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) + 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]));
    gsl_matrix_set(m, 2, 1, (-235.8 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) + 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]) + 1.143e4 * rho0 * (y[1] + y[0] + y[2]) * (0.0000134 + y[3]));
    gsl_matrix_set(m, 2, 2, (-235.8 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) - (471.7 * sqrt(rho0 * (y[1] + y[0] + y[2]))) / g0 + 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]));
    gsl_matrix_set(m, 2, 3, 1.143e4 * rho0 * y[1] * (y[1] + y[0] + y[2]));
    gsl_matrix_set(m, 2, 4, 0);

    gsl_matrix_set(m, 3, 0, (0.6175 * rho0 * y[2] * (0.0162 - 1. * y[3])) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))));
    gsl_matrix_set(m, 3, 1, (0.6175 * rho0 * y[2] * (0.0162 - 1. * y[3])) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))));
    gsl_matrix_set(m, 3, 2, (0.6175 * rho0 * y[2] * (0.0162 - 1. * y[3])) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) + (1.235 * sqrt(rho0 * (y[1] + y[0] + y[2])) * (0.0162 - 1. * y[3])) / g0);
    gsl_matrix_set(m, 3, 3, (-1.235 * y[2] * sqrt(rho0 * (y[1] + y[0] + y[2]))) / g0);
    gsl_matrix_set(m, 3, 4, 0);

    gsl_matrix_set(m, 4, 0, (0.5064 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))));
    gsl_matrix_set(m, 4, 1, (0.5064 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))));
    gsl_matrix_set(m, 4, 2, (0.5064 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) + (1.013 * sqrt(rho0 * (y[1] + y[0] + y[2]))) / g0);
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
    const double i0 = 0.9;
    const double a0 = 0.05;
    const double m0 = 0.05;
    const double z0 = 1e-4;
    const double s0 = 0.0;
    double Y0[] = {i0, a0, m0, z0, s0};

    /* Parameters */
    double params[] = {0.9 / 40.46, i0 + a0 + m0};

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

#ifndef EZ_SFR_H
#define EZ_SFR_H

/* Configuration */
#define MEMORY 1
#define IMF_PATH "./eta_data/Salpeter1955A/*" /* Initial mass function in use */
#define NUMEQU 5                               /* Number of equation in the model */

/* Error type constants */
#define GLOB_ERROR 0
#define GLS_ODE_ERROR 1
#define GSL_INT_ERROR 2

/* Unit conversion factors */
const double t_s = All.UnitTime_in_s / All.HubbleParam; /* T [internal_units] * t_s = T [s] */
const double t_Gyr = t_s / (SEC_PER_YEAR * 1e9);        /* T [internal_units] * t_Gyr = T [Gyr] */
const double factor[] = {
    8.411856872862986e-58,   /* factor[0] = f_ion = 1.0u"mp" -> u"Msun" */
    1.0094228247435582e-58}; /* factor[1] = f_diss = 0.4 * 0.15 * 2.0 * (1.0u"mp" -> u"Msun") */

/* ODE constants */
const double C1 = 0.809;    /* [Gyr Mₒ^(1/2) pc^(-3/2)] */
const double C2 = 3.024e-6; /* [Gyr Mₒ pc^(-3)] */
const double C3 = 8.745e-5; /* [Gyr Mₒ pc^(-3)] */
const double Zsun = 0.0134;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.145;
const double R = 0.1;

/* ODE integration constants */
const double abs_tol = 1.0e-8; /* Absolute tolerance */
const double rel_tol = 1.0e-8; /* Relative tolerance */

/* Density PDF from Burkhart (2018).
 *
 * Constant for the probability density function of gas density, using the
 * formulae from Burkhart et al. (2018) https://doi.org/10.3847/1538-4357/aad002.
 */
#define DIVISIONS 10
const double rho_pdf[DIVISIONS] = {
    0.12926709133405279, 0.19481409345869244,  0.21767875604365272,  0.18033768769448671,   0.11076747745448377,
    0.05043636613872607, 0.017021722247588605, 0.004256856424221968, 0.0007886571065396984, 0.0001127898050359872,
};
const double f_rho[DIVISIONS] = {
    0.0301973834223185, 0.0820849986238988, 0.22313016014842982, 0.6065306597126334, 1.6487212707001282,
    4.4816890703380645, 12.182493960703473, 33.11545195869231,   90.01713130052181,  244.69193226422038,
};

/* Structures */
struct Parameters
{
    double rho0;
    double g0;
    struct InterpFunc *interp_ion;
    struct InterpFunc *interp_diss;
};

struct InterpFunc
{
    gsl_spline *interp;
    gsl_interp_accel *acc;
    double x_min;
    double x_max;
};

// double **read_ftable(const char *filename, size_t *rows, size_t *cols);
// void check_error(const int err_type, const int err_code);
// double eval_interp(struct InterpFunc *interp, const double x);
// void get_eta(glob_t globbuf, const double age, const int ion, struct InterpFunc *interp_func);
// int sf_ode(double t, const double y[], double f[], void *params);
// int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
// void fractions(const double int_time, struct Parameters *params, double y[]);
// double stellar_fraction(const int index, const double int_time);

#endif /* #ifndef EZ_SFR_H */

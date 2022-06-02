#ifndef EZ_SFR_H
#define EZ_SFR_H

/* Initial mass function in use */
#define IMF_PATH "./eta_data/Salpeter1955B/*"

/* Error type constants */
#define GLOB_ERROR 0
#define GLS_ODE_ERROR 1
#define GSL_INT_ERROR 2

/* ODE constants */
#define const_1 0.809    /* [Gyr Mₒ^(1/2) pc^(-3/2)] */
#define const_2 3.024e-6 /* [Gyr Mₒ pc^(-3)] */
#define const_3 8.745e-5 /* [Gyr Mₒ pc^(-3)] */
#define Zsun 0.0134
#define Zeff (1e-3 * Zsun)
/*
 * Zsn and Rf were computed using the IMF from Salpeter et al. (1955)
 * in the mass range 0.1 M⊙ to 40 M⊙, as in Springel et al. (2003)
 */
#define Zsn 0.1157
#define Rf 0.0657

/* ODE integration constants */
#define ABS_TOL 1.0e-8 /* Absolute tolerance */
#define REL_TOL 0.0    /* Relative tolerance */

/* T [internal_units] * T_GYR = T [Gyr] */
#define T_GYR (All.UnitTime_in_s / All.HubbleParam / (SEC_PER_YEAR * 1e9))
/* RHO [internal_units] * RHO_CGS = RHO [g cm^(-3)] */
#define RHO_CGS (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv)
/* RHO [internal_units] * RHO_COSMO = RHO [Mₒ pc^(-3)] */
#define RHO_COSMO (RHO_CGS * PARSEC * PARSEC * PARSEC / SOLAR_MASS)
/* M [internal_units] * M_COSMO = M [Mₒ] */
#define M_COSMO (All.UnitMass_in_g / SOLAR_MASS)

/* 
 * Density PDF from Burkhart (2018) https://doi.org/10.3847/1538-4357/aad002 
 * Using 10 divisions in the range (-4, 6).
 */
#define DIVISIONS 10
static const double rho_pdf[] = {
    0.12926709133405279, 0.19481409345869244,  0.21767875604365272,  0.18033768769448671,   0.11076747745448377,
    0.05043636613872607, 0.017021722247588605, 0.004256856424221968, 0.0007886571065396984, 0.0001127898050359872,
};
static const double f_rho[] = {
    0.0301973834223185, 0.0820849986238988, 0.22313016014842982, 0.6065306597126334, 1.6487212707001282,
    4.4816890703380645, 12.182493960703473, 33.11545195869231,   90.01713130052181,  244.69193226422038,
};

/* Unit conversion constants */
/* T [internal_units] * T_GYR = T [Gyr] */
#define T_GYR (All.UnitTime_in_s / All.HubbleParam / (SEC_PER_YEAR * 1e9))
/* RHO [internal_units] * RHO_CGS = RHO [g cm^(-3)] */
#define RHO_CGS (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv)
/* RHO [internal_units] * RHO_COSMO = RHO [Mₒ pc^(-3)] */
#define RHO_COSMO (RHO_CGS * PARSEC * PARSEC * PARSEC / SOLAR_MASS)
/* M [internal_units] * M_COSMO = M [Mₒ] */
#define M_COSMO (All.UnitMass_in_g / SOLAR_MASS)

/* Structures */
struct ODEParameters
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
// int sf_ode(double t, const double y[], double f[], void *ode_params);
// int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *ode_params);
// void fractions(const double int_time, struct ODEParameters *ode_params, double y[]);
// double stellar_fraction(const int index, const double dt);
// void get_gas_fractions(const int index, double *fi, double *fa, double *fm);

#endif /* #ifndef EZ_SFR_H */

/*
 * Compile in Dirac with:
 * module swap gnu7/7.3.0 gnu10/10.3
 * module load gsl
 * gcc run_test.c -I/opt/ohpc/pub/libs/gnu10/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu10/gsl/2.6/lib -lgsl -lgslcblas -lm -o run_test
 */

#define _GNU_SOURCE

#include <glob.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>

#define GAMMA_MINUS1 (5.0 / 3) - 1
#define SOLAR_MASS 1.989e33
#define SEC_PER_YEAR 3.15576e7
#define PARSEC 3.085678e18

struct DD
{
    double Density;
};

struct ALL
{
    int ComovingIntegrationOn;
    double Time;
    double HubbleParam;
    double cf_a3inv;
    double MinEgySpec;
    double UnitDensity_in_cgs;
    double UnitTime_in_s;
};

struct SPHP
{
    double fHII; /* Ionized fraction */
    double fHI;  /* Atomic fraction */
    double fmol; /* Molecular fraction */
    double fstar;
    double fZ;
    double Entropy;
    double Ne;
    struct DD d;
};

void endrun(int a){};
double DMAX(double a, double b);
double AbundanceRatios(double a, double b, double *c, double *d, double *e);

int main(void)
{
    struct ALL All;

    All.ComovingIntegrationOn = 0;
    All.Time = 0.0;
    All.HubbleParam = 1.0;
    All.cf_a3inv = 1.0;
    All.MinEgySpec = 0.0;
    All.UnitDensity_in_cgs = 6.77e-22;
    All.UnitTime_in_s = 3.085678e16;

    struct SPHP SphP[1];

    SphP[0].fHII = 0.3;
    SphP[0].fHI = 0.3;
    SphP[0].fmol = 0.4;
    SphP[0].fstar = 0.0;
    SphP[0].fZ = 1e-4;
    SphP[0].Entropy = 0.0;
    SphP[0].Ne = 0.0;
    SphP[0].d.Density = 0.0022242332143521783;

#include "./ez_sfr.c"

    /***********************************************************************************************
     *
     * Test of the interpolation functions
     *
     **********************************************************************************************/

    printf("Test of the interpolation functions (eta_ion):\n\n");

    glob_t globbuf;
    int glob_err = glob(IMF_PATH, GLOB_ERR, NULL, &globbuf);

    struct InterpFunc interp_ion, interp_diss;

    interp_ion.acc = gsl_interp_accel_alloc();
    interp_ion.interp = gsl_spline_alloc(gsl_interp_linear, globbuf.gl_pathc);
    interp_ion.x_min = NAN;
    interp_ion.x_max = NAN;
    interp_diss.acc = gsl_interp_accel_alloc();
    interp_diss.interp = gsl_spline_alloc(gsl_interp_linear, globbuf.gl_pathc);
    interp_diss.x_min = NAN;
    interp_diss.x_max = NAN;

    get_eta(globbuf, 16.5, 0, &interp_ion);
    get_eta(globbuf, 16.5, 1, &interp_diss);

    printf(
        "eval_interp(&interp_ion, 0.0) = %g -- Should be %g\n",
        eval_interp(&interp_ion, 0.0),
        2498.936);
    printf(
        "eval_interp(&interp_ion, 0.0069) = %g -- Should be %g\n",
        eval_interp(&interp_ion, 0.0069),
        1325.183);
    printf(
        "eval_interp(&interp_ion, 0.004269) = %g -- Should be %g\n",
        eval_interp(&interp_ion, 0.004269),
        1486.551);
    printf(
        "eval_interp(&interp_ion, 0.69) = %g -- Should be %g\n\n",
        eval_interp(&interp_ion, 0.69),
        630.226);

    printf("Test of the interpolation functions (eta_diss):\n\n");

    printf(
        "eval_interp(&interp_diss, 0.0) = %g -- Should be %g\n",
        eval_interp(&interp_diss, 0.0),
        628.196);
    printf(
        "eval_interp(&interp_diss, 0.0069) = %g -- Should be %g\n",
        eval_interp(&interp_diss, 0.0069),
        451.739);
    printf(
        "eval_interp(&interp_diss, 0.004269) = %g -- Should be %g\n",
        eval_interp(&interp_diss, 0.004269),
        489.820);
    printf(
        "eval_interp(&interp_diss, 0.69) = %g -- Should be %g\n\n",
        eval_interp(&interp_diss, 0.69),
        291.946);

    /***********************************************************************************************
     *
     * Test of the ODEs system function
     *
     **********************************************************************************************/

    printf("Test of the ODEs system function:\n\n");

    double sol_ode[5] = {0};
    double sol_jac[25] = {0};
    double empty_vec[5] = {0};
    double ic[5] = {0.3, 0.3, 0.4, 1e-4, 0.0};
    double g0 = 1.0;
    struct ODEParameters parameters;
    double rho0 = 0.022243; // Msun * pc^-3

    parameters.rho0 = rho0;
    parameters.g0 = g0;
    parameters.interp_ion = &interp_ion;
    parameters.interp_diss = &interp_diss;

    sf_ode(0.0, ic, sol_ode, &parameters);

    printf("sol_ode[0] = %g -- Should be %g\n", sol_ode[0], -477.803);
    printf("sol_ode[1] = %g -- Should be %g\n", sol_ode[1], 524.114);
    printf("sol_ode[2] = %g -- Should be %g\n", sol_ode[2], -46.3801);
    printf("sol_ode[3] = %g -- Should be %g\n", sol_ode[3], 0.00055277);
    printf("sol_ode[4] = %g -- Should be %g\n\n", sol_ode[4], 0.0688853);

    /***********************************************************************************************
     *
     * Test of the jacobian function
     *
     **********************************************************************************************/

    printf("Test of the sol_jac function:\n\n");

    jacobian(0.0, ic, sol_jac, empty_vec, &parameters);

    printf("sol_jac[0] = %g -- Should be %g\n", sol_jac[0], -4321.52);
    printf("sol_jac[1] = %g -- Should be %g\n", sol_jac[1], 92.1219);
    printf("sol_jac[2] = %g -- Should be %g\n", sol_jac[2], 552.731);
    printf("sol_jac[3] = %g -- Should be %g\n", sol_jac[3], 0.0);
    printf("sol_jac[4] = %g -- Should be %g\n", sol_jac[4], 0.0);

    printf("sol_jac[5] = %g -- Should be %g\n", sol_jac[5], 4344.68);
    printf("sol_jac[6] = %g -- Should be %g\n", sol_jac[6], -68.9995);
    printf("sol_jac[7] = %g -- Should be %g\n", sol_jac[7], -413.781);
    printf("sol_jac[8] = %g -- Should be %g\n", sol_jac[8], -76.3035);
    printf("sol_jac[9] = %g -- Should be %g\n", sol_jac[9], 0.0);

    printf("sol_jac[10] = %g -- Should be %g\n", sol_jac[10], -23.1857);
    printf("sol_jac[11] = %g -- Should be %g\n", sol_jac[11], -23.1569);
    printf("sol_jac[12] = %g -- Should be %g\n", sol_jac[12], -139.158);
    printf("sol_jac[13] = %g -- Should be %g\n", sol_jac[13], 76.3035);
    printf("sol_jac[14] = %g -- Should be %g\n", sol_jac[14], 0.0);

    printf("sol_jac[15] = %g -- Should be %g\n", sol_jac[15], 0.000276385);
    printf("sol_jac[16] = %g -- Should be %g\n", sol_jac[16], 0.000276385);
    printf("sol_jac[17] = %g -- Should be %g\n", sol_jac[17], 0.00165831);
    printf("sol_jac[18] = %g -- Should be %g\n", sol_jac[18], -0.073727);
    printf("sol_jac[19] = %g -- Should be %g\n", sol_jac[19], 0.0);

    printf("sol_jac[20] = %g -- Should be %g\n", sol_jac[20], 0.0344426);
    printf("sol_jac[21] = %g -- Should be %g\n", sol_jac[21], 0.0344426);
    printf("sol_jac[22] = %g -- Should be %g\n", sol_jac[22], 0.206656);
    printf("sol_jac[23] = %g -- Should be %g\n", sol_jac[23], 0.0);
    printf("sol_jac[24] = %g -- Should be %g\n\n", sol_jac[24], 0.0);

    /***********************************************************************************************
     *
     * Test of the function stellar_fraction()
     *
     **********************************************************************************************/

    printf("Test of the function stellar_fraction():\n\n");

    double sf = stellar_fraction(0, 1.0227120263358653);

    printf("star fraction = %g -- Should be %g\n", sf, 0.000586533712524996);

    return 0;
}

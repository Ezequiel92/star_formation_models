#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/*******************************************************************************************
* Model 005
*
* Alternative names: 'Modified basic model with a non-constant number density'
*******************************************************************************************/

static const char NAME[] = "Model 005";     // Model name
static const int SIZE = 5;                  // Size of the system of ODEs

/* File initialization function */

FILE *file_ini(const char *name)
{
    FILE *file;
    file = fopen(name, "a");
    fprintf(file, "t \t i(t) \t a(t) \t m(t) \t z(t) \t s(t) \n");
    return file;
}

/* Constants. */

const double K1 = 0.001656;			// [pc^4 Mₒ^(-2) cm^(-3)]
const double K2 = 0.00984;			// [pc^4 Mₒ^(-2) cm^(-3)]
const double tau_S = 2.6;       	// [Gyr]
const double Rh = 1.9;          	// [cm^3 Gyr^-1]
const double SigmaNuB = 8158.0; 	// [cm^3 Gyr^-1]
const double Zsun = 0.0134;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.09;
const double eta_ion = 955.29;
const double eta_diss = 380.93;
const double R = 0.18;			

/* System of equations */

/*
	Ionized gas fraction:       i(t) / g -> y[0]
    Atomic gas fraction:        a(t) / g -> y[1]
    Molecular gas fraction:     m(t) / g -> y[2]
    Metal fraction:             z(t) / g -> y[3]
    Stellar fraction:           s(t) / g -> y[4]	

    where g = i(t) + a(t) + m(t)
     
    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
*/
static int ode_system(double t, const double y[], double f[], void *params)
{
	(void)(t);
	
	/* Parameters */
	double g = *(double *)params;
	double Z = *((double *)params + 1);
	
	/* Auxiliary equations */
	double ne = K1 * (1 + y[3]) * g * g;
	double nh = K2 * (1 + y[3]) * g * g;
    double tau_R = 1 / (ne * SigmaNuB);
    double tau_C = Zsun / (2 * nh * Rh * (Z + Zeff));
    double recombination = y[0] / tau_R;
    double cloud_formation = y[1] / tau_C;
    double psi = y[2] / tau_S;
	
	/* ODE system */
	f[0] = -recombination + eta_ion * psi;
    f[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
    f[2] = cloud_formation - (1 + eta_diss) * psi;
    f[3] = (Zsn * R - Z) * psi;
	f[4] = psi;

	return GSL_SUCCESS;
};

/* Jacobian */

int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	(void)(t);

	double g = *(double *)params;
	double Z = *((double *)params + 1);

	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 5, 5);
	gsl_matrix * m = &dfdy_mat.matrix;

	gsl_matrix_set(m, 0, 0, -13.509648*pow(g,2)*(1. + y[3]));
	gsl_matrix_set(m, 0, 1, 0);
	gsl_matrix_set(m, 0, 2, 367.4192308);
	gsl_matrix_set(m, 0, 3, -13.509648*pow(g,2)*y[0]);
	gsl_matrix_set(m, 0, 4, 0);

	gsl_matrix_set(m, 1, 0, 13.509648*pow(g,2)*(1. + y[3]));
	gsl_matrix_set(m, 1, 1, -2.790447761*pow(g,2)*(0.0000134 + Z)*(1. + y[3]));
	gsl_matrix_set(m, 1, 2, -220.9076923);
	gsl_matrix_set(m, 1, 3, -2.790447761*pow(g,2)*(0.0000134 + Z)*y[1] + 13.509648*pow(g,2)*y[0]);
	gsl_matrix_set(m, 1, 4, 0);

	gsl_matrix_set(m, 2, 0, 0);
	gsl_matrix_set(m, 2, 1, 2.790447761*pow(g,2)*(0.0000134 + Z)*(1. + y[3]));
	gsl_matrix_set(m, 2, 2, -146.8961538);
	gsl_matrix_set(m, 2, 3, 2.790447761*pow(g,2)*(0.0000134 + Z)*y[1]);
	gsl_matrix_set(m, 2, 4, 0);

	gsl_matrix_set(m, 3, 0, 0);
	gsl_matrix_set(m, 3, 1, 0);
	gsl_matrix_set(m, 3, 2, 0.3846153846*(0.0162 - 1.*Z));
	gsl_matrix_set(m, 3, 3, 0);
	gsl_matrix_set(m, 3, 4, 0);

	gsl_matrix_set(m, 4, 0, 0);
	gsl_matrix_set(m, 4, 1, 0);
	gsl_matrix_set(m, 4, 2, 0.3846153846);
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

double t_start = 0.0;    // [Gyr]
double t_end = 1.0;      // [Gyr]

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
	double params[] = {1e-3 * 14.8285, z0};       // {g [Mₒ pc^(-2)], Z}

	/* Integration */
	gsl_odeiv2_system sys = {ode_system, jacobian, SIZE, params};
	gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
        &sys, 
        gsl_odeiv2_step_msbdf, 
        dt, 
        abs_tol, 
        rel_tol
    );

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
            t_start, Y0[0], Y0[1], Y0[2], Y0[3], Y0[4]
        );
	}

	gsl_odeiv2_driver_free(driver);
	fclose(file);

    printf("Stellar fraction after %gGyr (msbdf): %.3g", t_end, Y0[4]);
}


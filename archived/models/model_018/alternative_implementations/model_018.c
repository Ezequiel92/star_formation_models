#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/*******************************************************************************************
* Model 018
*
* Alternative names: 'Fully coupled model (constant τS)', 'Model E'
*******************************************************************************************/

static const char NAME[] = "Model 018";     // Model name
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

const double tau_S = 2.6;       	// [Gyr]
const double C1 = 3.03e-6; // [Gyr]
const double C2 = 8.74e-5; // [Gyr]
const double P = 40.4;	   // [cm^(-3) Mₒ^(-1) pc^3]
const double Zsun = 0.0134;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.09;
const double eta_ion = 955.29;
const double eta_diss = 380.93;
const double R = 0.18;			

/* System of equations */

/*
	Ionized gas fraction:       i(t) / tot -> y[0]
    Atomic gas fraction:        a(t) / tot -> y[1]
    Molecular gas fraction:     m(t) / tot -> y[2]
    Metal fraction:             z(t) / tot -> y[3]
    Stellar fraction:           s(t) / tot -> y[4]	

    where tot = i(t) + a(t) + m(t) + z(t) + s(t)
     
    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
*/
static int ode_system(double t, const double y[], double f[], void *params)
{
	(void)(t);
	
	/* Parameters */
	double tot0 = *(double *)params;
	
	/* Auxiliary equations */
	double g = y[0] + y[1] + y[2] + y[3];
	double psi = y[2] / tau_S;
	double ig = y[0] / g;
	double ag = y[1] / g;
	double mg = y[2] / g;
	double zg = y[3] / g;
    double tau_R = C1 / (y[0] * tot0);
	double tau_C = C2 / ((y[1] + y[2]) * tot0 * (zg + Zeff));
    double recombination = y[0] / tau_R;
	double cloud_formation = y[1] / tau_C;
	
	/* ODE system */
	f[0] = -recombination + (eta_ion + (1 - Zsn) * R - ig) * psi;
    f[1] = -cloud_formation + recombination + (eta_diss - eta_ion - ag) * psi;
    f[2] = cloud_formation - (eta_diss + mg) * psi;
    f[3] = (Zsn * R - zg) * psi;
	f[4] = (1 - R) * psi;

	return GSL_SUCCESS;
};

/* Jacobian */

int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	(void)(t);

	double tot0 = *(double *)params;

	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 5, 5);
	gsl_matrix * m = &dfdy_mat.matrix;

	gsl_matrix_set(m, 0, 0, -660066.0066*tot0*y[0] + 0.3846153846*y[2]*(y[0]/pow(y[1] + y[0] + y[2] + y[3],2) - 1./(y[1] + y[0] + y[2] + y[3])));
	gsl_matrix_set(m, 0, 1, (0.3846153846*y[0]*y[2])/pow(y[1] + y[0] + y[2] + y[3],2));
	gsl_matrix_set(m, 0, 2, (0.3846153846*y[0]*y[2])/pow(y[1] + y[0] + y[2] + y[3],2) + 0.3846153846*(955.4538 - (1.*y[0])/(y[1] + y[0] + y[2] + y[3])));
	gsl_matrix_set(m, 0, 3, (0.3846153846*y[0]*y[2])/pow(y[1] + y[0] + y[2] + y[3],2));
	gsl_matrix_set(m, 0, 4, 0);

	gsl_matrix_set(m, 1, 0, 660066.0066*tot0*y[0] + (0.3846153846*y[1]*y[2])/pow(y[1] + y[0] + y[2] + y[3],2) + (11441.6476*tot0*y[1]*(y[1] + y[2])*y[3])/pow(y[1] + y[0] + y[2] + y[3],2));
	gsl_matrix_set(m, 1, 1, (11441.6476*tot0*y[1]*(y[1] + y[2])*y[3])/pow(y[1] + y[0] + y[2] + y[3],2) + 0.3846153846*y[2]*(y[1]/pow(y[1] + y[0] + y[2] + y[3],2) - 1./(y[1] + y[0] + y[2] + y[3])) - 11441.6476*tot0*y[1]*(0.0000134 + y[3]/(y[1] + y[0] + y[2] + y[3])) - 11441.6476*tot0*(y[1] + y[2])*(0.0000134 + y[3]/(y[1] + y[0] + y[2] + y[3])));
	gsl_matrix_set(m, 1, 2, (0.3846153846*y[1]*y[2])/pow(y[1] + y[0] + y[2] + y[3],2) + (11441.6476*tot0*y[1]*(y[1] + y[2])*y[3])/pow(y[1] + y[0] + y[2] + y[3],2) + 0.3846153846*(-574.36 - (1.*y[1])/(y[1] + y[0] + y[2] + y[3])) - 11441.6476*tot0*y[1]*(0.0000134 + y[3]/(y[1] + y[0] + y[2] + y[3])));
	gsl_matrix_set(m, 1, 3, (0.3846153846*y[1]*y[2])/pow(y[1] + y[0] + y[2] + y[3],2) - 11441.6476*tot0*y[1]*(y[1] + y[2])*((-1.*y[3])/pow(y[1] + y[0] + y[2] + y[3],2) + 1/(y[1] + y[0] + y[2] + y[3])));
	gsl_matrix_set(m, 1, 4, 0);

	gsl_matrix_set(m, 2, 0, (0.3846153846*pow(y[2],2))/pow(y[1] + y[0] + y[2] + y[3],2) - (11441.6476*tot0*y[1]*(y[1] + y[2])*y[3])/pow(y[1] + y[0] + y[2] + y[3],2));
	gsl_matrix_set(m, 2, 1, (0.3846153846*pow(y[2],2))/pow(y[1] + y[0] + y[2] + y[3],2) - (11441.6476*tot0*y[1]*(y[1] + y[2])*y[3])/pow(y[1] + y[0] + y[2] + y[3],2) + 11441.6476*tot0*y[1]*(0.0000134 + y[3]/(y[1] + y[0] + y[2] + y[3])) + 11441.6476*tot0*(y[1] + y[2])*(0.0000134 + y[3]/(y[1] + y[0] + y[2] + y[3])));
	gsl_matrix_set(m, 2, 2, (-11441.6476*tot0*y[1]*(y[1] + y[2])*y[3])/pow(y[1] + y[0] + y[2] + y[3],2) - 0.3846153846*y[2]*((-1.*y[2])/pow(y[1] + y[0] + y[2] + y[3],2) + 1/(y[1] + y[0] + y[2] + y[3])) - 0.3846153846*(380.93 + y[2]/(y[1] + y[0] + y[2] + y[3])) + 11441.6476*tot0*y[1]*(0.0000134 + y[3]/(y[1] + y[0] + y[2] + y[3])));
	gsl_matrix_set(m, 2, 3, (0.3846153846*pow(y[2],2))/pow(y[1] + y[0] + y[2] + y[3],2) + 11441.6476*tot0*y[1]*(y[1] + y[2])*((-1.*y[3])/pow(y[1] + y[0] + y[2] + y[3],2) + 1/(y[1] + y[0] + y[2] + y[3])));
	gsl_matrix_set(m, 2, 4, 0);

	gsl_matrix_set(m, 3, 0, (0.3846153846*y[2]*y[3])/pow(y[1] + y[0] + y[2] + y[3],2));
	gsl_matrix_set(m, 3, 1, (0.3846153846*y[2]*y[3])/pow(y[1] + y[0] + y[2] + y[3],2));
	gsl_matrix_set(m, 3, 2, (0.3846153846*y[2]*y[3])/pow(y[1] + y[0] + y[2] + y[3],2) + 0.3846153846*(0.0162 - (1.*y[3])/(y[1] + y[0] + y[2] + y[3])));
	gsl_matrix_set(m, 3, 3, 0.3846153846*y[2]*(y[3]/pow(y[1] + y[0] + y[2] + y[3],2) - 1./(y[1] + y[0] + y[2] + y[3])));
	gsl_matrix_set(m, 3, 4, 0);

	gsl_matrix_set(m, 4, 0, 0);
	gsl_matrix_set(m, 4, 1, 0);
	gsl_matrix_set(m, 4, 2, 0.3153846154);
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
	double params[] = {1e-3 / P};	// {tot0 [Mₒ pc^(-3)] = n / P}

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


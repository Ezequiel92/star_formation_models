#include <cstdio>
#include <cmath>
#include <cstring>
#include <boost/numeric/odeint.hpp>

typedef double Double;
typedef boost::numeric::ublas::vector<Double> vector_type;
typedef boost::numeric::ublas::matrix<Double> matrix_type;

using namespace std;
using namespace boost::numeric::odeint;

/* Parameters used within the models to be varied in the main loop. */

Double n; // [1/cm^3]
Double g; // [Mₒ pc^(-2)]
Double Z; // [Mₒ pc^(-2)]

const Double mass_factor = 14.8285; // (proton_mass / cm^3) * 0.6 kpc = 14.8285 Mₒ / pc^2 
									// where 0.6 kpc is in order of magnitude the height of the galaxy

/***************************************************************	
	Ranges of parameters and initial conditions.
	x_min = 10^min_x, x_max = 10^max_x, dx = 10^delta_x. 
****************************************************************/

const Double min_a = -4;
const Double max_a = 0;
const Double delta_a = 1.5;

const Double min_m = -4.;
const Double max_m = 0;
const Double delta_m = 1.5;

const Double min_z = -4.;
const Double max_z = -1;
const Double delta_z = 1.5;

const Double min_n = -1;
const Double max_n = 3;
const Double delta_n = 1.5;

/* Number of values for each parameter to be swept in the loop. */

const int num_a = (max_a - min_a) / delta_a + 1;
const int num_m = (max_m - min_m) / delta_m + 1;
const int num_z = (max_z - min_z) / delta_z + 1;
const int num_n = (max_n - min_n) / delta_n + 1;

/* Models. */

//#include "models/basic_model.hpp"
//#include "models/full_model_constant_tau_no_IO.hpp"
//#include "models/full_model_varying_tau_no_IO.hpp"
//#include "models/simplified_model.hpp"
//#include "models/full_model_kennicutt_no_IO.hpp"
#include "models/ascasibar_model_constant_tau.hpp"

int main()
{
	const Double Tstart = 0.0;				  // Integration initial time [Gyr].
	const Double Tend = 10.0;	  			  // Integration final time [Gyr].
	const Double dt = (Tend - Tstart) * 1e-5; // Time step.
	Double abs_err = 1.0e-8;				  // Absolute error.
	Double rel_err = 1.0e-8;				  // Relative error.

	FILE *file;
	char filepath[90];

	/* Integration with Rosenbrock4. */
	auto stepper_rosenbrock = make_dense_output(abs_err, rel_err, rosenbrock4<Double>());

	if (strcmp(name, "basic_model") == 0)
	{
		for (size_t i_n = 0; i_n < num_n; ++i_n)
		{
			for (size_t i_a = 0; i_a < num_a; ++i_a)
			{
				for (size_t i_m = 0; i_m < num_m; ++i_m)
				{
					for (size_t i_z = 0; i_z < num_z; ++i_z)
					{
						n = pow(10., (min_n + i_n * delta_n));		 // [1/cm^3]
						Y0[1] = pow(10., (min_a + i_a * delta_a));
						Y0[2] = pow(10., (min_m + i_m * delta_m));
						Y0[3] = pow(10., (min_z + i_z * delta_z));
						Z = Y0[3];

						if (Y0[1] + Y0[2] <= 1)
						{
							Y0[0] = 1 - Y0[1] - Y0[2];
							sprintf(filepath, "results/data/%s_%d_%d_%d_%d.dat",
									name, i_n, i_a, i_m, i_z);
							file = fopen(filepath, "a");
							fprintf(file, "n = %.3g [cm^-3]  if0 = %.3g  af0 = %.3g  mf0 = %.3g  zf0 = %.3g\n",
									n, Y0[0], Y0[1], Y0[2], Y0[3]);
							integrate_const(stepper_rosenbrock, make_pair(rhs(), jacobi()), Y0, Tstart, Tend, dt, file_observer(file));
							fclose(file);
						}
					}
				}
			}
		}
	}
	else if (strcmp(name, "simplified_model") == 0)
	{
		for (size_t i_n = 0; i_n < num_n; ++i_n)
		{
			for (size_t i_a = 0; i_a < num_a; ++i_a)
			{
				for (size_t i_m = 0; i_m < num_m; ++i_m)
				{
					for (size_t i_z = 0; i_z < num_z; ++i_z)
					{
						g = mass_factor * pow(10., (min_n + i_n * delta_n));	// [Mₒ / pc^2]
						Y0[1] = pow(10., (min_a + i_a * delta_a));
						Y0[2] = pow(10., (min_m + i_m * delta_m));
						Y0[3] = 0.0;
						Y0[4] = pow(10., (min_z + i_z * delta_z));

						if (Y0[1] + Y0[2] <= 1)
						{
							Y0[0] = 1 - Y0[1] - Y0[2];
							sprintf(filepath, "results/data/%s_%d_%d_%d_%d.dat",
									name, i_n, i_a, i_m, i_z);
							file = fopen(filepath, "a");
							fprintf(file, "g = %.3g [Mₒ / pc^2]  if0 = %.3g  af0 = %.3g  mf0 = %.3g  zf0 = %.3g\n",
									g, Y0[0], Y0[1], Y0[2], Y0[3], Y0[4]);
							integrate_const(stepper_rosenbrock, make_pair(rhs(), jacobi()), Y0, Tstart, Tend, dt, file_observer(file));
							fclose(file);
						}
					}
				}
			}
		}
	}
	else
	{
		for (size_t i_n = 0; i_n < num_n; ++i_n)
		{
			for (size_t i_a = 0; i_a < num_a; ++i_a)
			{
				for (size_t i_m = 0; i_m < num_m; ++i_m)
				{
					for (size_t i_z = 0; i_z < num_z; ++i_z)
					{
						g = mass_factor * pow(10., (min_n + i_n * delta_n));	// [Mₒ / pc^2]
						Y0[1] = pow(10., (min_a + i_a * delta_a));
						Y0[2] = pow(10., (min_m + i_m * delta_m));
						Y0[3] = 0.0;
						Y0[4] = pow(10., (min_z + i_z * delta_z));

						if (Y0[1] + Y0[2] <= 1)
						{
							Y0[0] = g * (1 - Y0[1] - Y0[2]);	// [Mₒ / pc^2]
							Y0[1] *= g;							// [Mₒ / pc^2]
							Y0[2] *= g;							// [Mₒ / pc^2]
							Y0[3] *= g;							// [Mₒ / pc^2]
							Y0[4] *= g;							// [Mₒ / pc^2]
							sprintf(filepath, "results/data/%s_%d_%d_%d_%d.dat",
									name, i_n, i_a, i_m, i_z);
							file = fopen(filepath, "a");
							fprintf(file, "g = %.3g  i0 = %.3g  a0 = %.3g  m0 = %.3g  z0 = %.3g  --  All in [Mₒ / pc^2]\n",
									g, Y0[0], Y0[1], Y0[2], Y0[3], Y0[4]);
							integrate_const(stepper_rosenbrock, make_pair(rhs(), jacobi()), Y0, Tstart, Tend, dt, file_observer(file));
							fclose(file);
						}
					}
				}
			}
		}
	}
}

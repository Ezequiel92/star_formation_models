#include <cstdio>
#include <cmath>
#include <chrono>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <boost/numeric/odeint.hpp>

typedef double Double;
typedef boost::numeric::ublas::vector<Double> vector_type;
typedef boost::numeric::ublas::matrix<Double> matrix_type;

using namespace std;
using namespace std::chrono;
using namespace boost::numeric::odeint;

Double Mean(Double *numbers, size_t length)
{
	/* Mean value of 'numbers'. */

	return (Double)accumulate(numbers, numbers + length, 0.0) / length;
}

Double StandardDeviation(Double *numbers, Double mean, size_t length)
{
	/* Standard deviation of 'numbers' with mean value 'mean'. */

	Double sDeviation = 0.0;
	for (size_t i = 0; i < length; ++i)
	{
		sDeviation += pow(numbers[i] - mean, 2);
	}

	return sqrt(sDeviation / (length - 1));
}

/* Models. */

//#include "models/basic_model.hpp"
//#include "models/full_model_constant_tau.hpp"
//#include "models/full_model_constant_tau_no_IO.hpp"
//#include "models/full_model_varying_tau_no_IO.hpp"
//#include "models/simplified_model.hpp"
//#include "models/full_model_kennicutt_no_IO.hpp"
//#include "models/full_model_varying_tau.hpp"
//#include "models/ascasibar_model_constant_tau.hpp"
//#include "models/modified_basic_model.hpp"
//#include "models/modified_basic_model_coupled_metallicity.hpp"
//#include "models/modified_basic_model_mass_recycling.hpp"
//#include "models/modified_basic_model_non-constant_number_density.hpp"
//#include "models/modified_basic_model_three_changes.hpp"
#include "models/self_consistent_ascasibar_model_constant_tau.hpp"

int main()
{
	const Double dt = (Tend - Tstart) * 1e-5; // Time step.
	Double abs_err = 1.0e-6;				  // Absolute error.
	Double rel_err = 1.0e-8;				  // Relative error.

	size_t n_steps; // Number of steps made by the numerical integration method.
	size_t n = 100; // Repetitions of the time measurement to get good statistics.

	/* Variables to save timing results. */
	Double times[n]{};
	steady_clock::time_point begin, end;
	Double time_span;

	/* File to store results */
	FILE *file;
	file = fopen("results/time_comparison_models.dat", "a");

	if (strcmp(name, "Basic model") == 0)
	{
		fprintf(file, "<Model>\t\t<Total duration>\t\t<Time per step>\t\t<Number of steps> \n\n");

		/* Integration with Rosenbrock4. */
		auto stepper_rosenbrock = make_dense_output(abs_err, rel_err, rosenbrock4<Double>());
		for (size_t i = 0; i < n; ++i)
		{
			Y0[0] = i0;
			Y0[1] = a0;
			Y0[2] = m0;
			Y0[3] = z0;
			begin = steady_clock::now();
			n_steps = integrate_const(stepper_rosenbrock, make_pair(rhs(), jacobi_const()), Y0, Tstart, Tend, dt);
			end = steady_clock::now();
			time_span = duration_cast<duration<Double, micro>>(end - begin).count();
			times[i] = time_span;
		}

		/* Mean and standard deviation of the time taken by one full integration in milliseconds. */
		Double mean_time = Mean(times, n) / 1e3;
		Double std_time = StandardDeviation(times, mean_time * 1e3, n) / 1e3;

		/* Mean and standard deviation of time taken by one integration step in milliseconds. */
		Double mean_time_per_step = mean_time / (Double)n_steps;
		Double std_time_per_step = mean_time_per_step * (std_time / mean_time);

		/* Save results in file. */
		fprintf(file, "%s:\t\t(%.2g ± %.1g) ms\t\t(%.2g ± %.1g) ms/step\t\t%zu steps \n\n",
				name, mean_time, std_time, mean_time_per_step, std_time_per_step, n_steps);
		fclose(file);
	}
	else if (strcmp(name, "Modified basic model") == 0 || strcmp(name, "Modified basic model with mass recycling") == 0)
	{
		/* Integration with Rosenbrock4. */
		auto stepper_rosenbrock = make_dense_output(abs_err, rel_err, rosenbrock4<Double>());
		for (size_t i = 0; i < n; ++i)
		{
			Y0[0] = i0;
			Y0[1] = a0;
			Y0[2] = m0;
			Y0[3] = s0;
			Y0[4] = z0;
			begin = steady_clock::now();
			n_steps = integrate_const(stepper_rosenbrock, make_pair(rhs(), jacobi_const()), Y0, Tstart, Tend, dt);
			end = steady_clock::now();
			time_span = duration_cast<duration<Double, micro>>(end - begin).count();
			times[i] = time_span;
		}

		/* Mean and standard deviation of the time taken by one full integration in milliseconds. */
		Double mean_time = Mean(times, n) / 1e3;
		Double std_time = StandardDeviation(times, mean_time * 1e3, n) / 1e3;

		/* Mean and standard deviation of time taken by one integration step in milliseconds. */
		Double mean_time_per_step = mean_time / (Double)n_steps;
		Double std_time_per_step = mean_time_per_step * (std_time / mean_time);

		/* Save results in file. */
		fprintf(file, "%s:\t\t(%.2g ± %.1g) ms\t\t(%.2g ± %.1g) ms/step\t\t%zu steps \n\n",
				name, mean_time, std_time, mean_time_per_step, std_time_per_step, n_steps);
		fclose(file);
	}
	else
	{
		/* Integration with Rosenbrock4. */
		auto stepper_rosenbrock = make_dense_output(abs_err, rel_err, rosenbrock4<Double>());
		for (size_t i = 0; i < n; ++i)
		{
			Y0[0] = i0;
			Y0[1] = a0;
			Y0[2] = m0;
			Y0[3] = s0;
			Y0[4] = z0;
			begin = steady_clock::now();
			n_steps = integrate_const(stepper_rosenbrock, make_pair(rhs(), jacobi_var()), Y0, Tstart, Tend, dt);
			end = steady_clock::now();
			time_span = duration_cast<duration<Double, micro>>(end - begin).count();
			times[i] = time_span;
		}

		/* Mean and standard deviation of the time taken by one full integration in milliseconds. */
		Double mean_time = Mean(times, n) / 1e3;
		Double std_time = StandardDeviation(times, mean_time * 1e3, n) / 1e3;

		/* Mean and standard deviation of time taken by one integration step in milliseconds. */
		Double mean_time_per_step = mean_time / (Double)n_steps;
		Double std_time_per_step = mean_time_per_step * (std_time / mean_time);

		/* Save results in file. */
		fprintf(file, "%s:\t\t(%.2g ± %.1g) ms\t\t(%.2g ± %.1g) ms/step\t\t%zu steps \n\n",
				name, mean_time, std_time, mean_time_per_step, std_time_per_step, n_steps);
		fclose(file);
	}
}

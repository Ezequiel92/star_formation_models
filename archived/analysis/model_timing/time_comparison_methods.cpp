 see#include <cstdio>
#include <cmath>
#include <chrono>
#include <boost/numeric/odeint.hpp>

typedef double Double;
typedef boost::numeric::ublas::vector<Double> vector_type;
typedef boost::numeric::ublas::matrix<Double> matrix_type;

using namespace std;
using namespace std::chrono;
using namespace boost::numeric::odeint;

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

	size_t n = 50; // Number of iterations in for loop.

	/* Number of steps each methods requires. */
	size_t n_steps_rosenbrock_CJ;
	size_t n_steps_rosenbrock_VJ;
	size_t n_steps_rkdopri;

	/* File to store results */
	FILE *file;
	file = fopen("results/time_comparison_methods.dat", "a");

	if (strcmp(name, "Basic model") == 0)
	{
		/* Integration with Rosenbrock4. */
		auto stepper_rosenbrock = make_dense_output(abs_err, rel_err, rosenbrock4<Double>());
		auto begin = steady_clock::now();
		for (size_t i = 0; i < n; ++i)
		{
			Y0[0] = i0;
			Y0[1] = a0;
			Y0[2] = m0;
			Y0[3] = z0;
			n_steps_rosenbrock_CJ = integrate_const(stepper_rosenbrock, make_pair(rhs(), jacobi_const()), Y0, Tstart, Tend, dt);
		}
		auto end = steady_clock::now();

		// Total time taken for n iterations in nanoseconds:
		auto duration_rosenbrock_CJ = duration_cast<duration<Double, nano>>(end - begin).count();
		// Time taken by one iteration in seconds:
		duration_rosenbrock_CJ = duration_rosenbrock_CJ / (1e9 * n);
		// Time taken by one step in milliseconds:
		Double durationPS_rosenbrock_CJ = (duration_rosenbrock_CJ / (Double)n_steps_rosenbrock_CJ) * 1e3;

		/* Integration with Runge-Kutta-Dormand–Prince. */
		auto stepper_rkdopri = make_dense_output(abs_err, rel_err, runge_kutta_dopri5<vector_type>());
		begin = steady_clock::now();
		for (size_t i = 0; i < n; ++i)
		{
			Y0[0] = i0;
			Y0[1] = a0;
			Y0[2] = m0;
			Y0[3] = z0;
			n_steps_rkdopri = integrate_const(stepper_rkdopri, rhs(), Y0, Tstart, Tend, dt);
		}
		end = steady_clock::now();

		// Total time taken for n iterations in nanoseconds:
		auto duration_rkdopri = duration_cast<duration<Double, nano>>(end - begin).count();
		// Time taken by one iteration in seconds:
		duration_rkdopri = duration_rkdopri / (1e9 * n);
		// Time taken by one step in milliseconds:
		Double durationPS_rkdopri = (duration_rkdopri / (Double)n_steps_rkdopri) * 1e3;

		/* Save results in file. */

		Double duration[] = {duration_rosenbrock_CJ, duration_rkdopri};
		auto min_dur = *min_element(duration, duration + 2);

		size_t steps[] = {n_steps_rosenbrock_CJ, n_steps_rkdopri};
		auto min_steps = (Double)*min_element(steps, steps + 2);

		Double time_step[] = {durationPS_rosenbrock_CJ, durationPS_rkdopri};
		auto min_time_step = *min_element(time_step, time_step + 2);

		fprintf(file, "%s:\t\t<Total duration>\t\t<Time per step>\t\t<Number of steps> \n\n", name);
		fprintf(file, "\t\tRosenbrock (constant Jacobian):\t\t%.2g s (x%.2g)\t\t%.2g ms/step (x%.2g)\t\t%zu (x%.2g) \n",
				duration_rosenbrock_CJ, duration_rosenbrock_CJ / min_dur,
				durationPS_rosenbrock_CJ, durationPS_rosenbrock_CJ / min_time_step,
				n_steps_rosenbrock_CJ, (Double)n_steps_rosenbrock_CJ / min_steps);
		fprintf(file, "\t\tRKDOPRI:                       \t\t%.2g s (x%.2g)\t\t%.2g ms/step (x%.2g)\t\t%zu (x%.2g) \n\n",
				duration_rkdopri, duration_rkdopri / min_dur,
				durationPS_rkdopri, durationPS_rkdopri / min_time_step,
				n_steps_rkdopri, (Double)n_steps_rkdopri / min_steps);
		fclose(file);
	}
	else if (strcmp(name, "Modified basic model") == 0 || strcmp(name, "Modified basic model with mass recycling") == 0)
	{
		/* Integration with Rosenbrock4. */
		auto stepper_rosenbrock = make_dense_output(abs_err, rel_err, rosenbrock4<Double>());
		auto begin = steady_clock::now();
		for (size_t i = 0; i < n; ++i)
		{
			Y0[0] = i0;
			Y0[1] = a0;
			Y0[2] = m0;
			Y0[3] = s0;
			Y0[4] = z0;
			n_steps_rosenbrock_CJ = integrate_const(stepper_rosenbrock, make_pair(rhs(), jacobi_const()), Y0, Tstart, Tend, dt);
		}
		auto end = steady_clock::now();

		// Total time taken for n iterations in nanoseconds:
		auto duration_rosenbrock_CJ = duration_cast<duration<Double, nano>>(end - begin).count();
		// Time taken by one iteration in seconds:
		duration_rosenbrock_CJ = duration_rosenbrock_CJ / (1e9 * n);
		// Time taken by one step in milliseconds:
		Double durationPS_rosenbrock_CJ = (duration_rosenbrock_CJ / (Double)n_steps_rosenbrock_CJ) * 1e3;

		/* Integration with Runge-Kutta-Dormand–Prince. */
		auto stepper_rkdopri = make_dense_output(abs_err, rel_err, runge_kutta_dopri5<vector_type>());
		begin = steady_clock::now();
		for (size_t i = 0; i < n; ++i)
		{
			Y0[0] = i0;
			Y0[1] = a0;
			Y0[2] = m0;
			Y0[3] = s0;
			Y0[4] = z0;
			n_steps_rkdopri = integrate_const(stepper_rkdopri, rhs(), Y0, Tstart, Tend, dt);
		}
		end = steady_clock::now();

		// Total time taken for n iterations in nanoseconds:
		auto duration_rkdopri = duration_cast<duration<Double, nano>>(end - begin).count();
		// Time taken by one iteration in seconds:
		duration_rkdopri = duration_rkdopri / (1e9 * n);
		// Time taken by one step in milliseconds:
		Double durationPS_rkdopri = (duration_rkdopri / (Double)n_steps_rkdopri) * 1e3;

		/* Save results in file. */

		Double duration[] = {duration_rosenbrock_CJ, duration_rkdopri};
		auto min_dur = *min_element(duration, duration + 2);

		size_t steps[] = {n_steps_rosenbrock_CJ, n_steps_rkdopri};
		auto min_steps = (Double)*min_element(steps, steps + 2);

		Double time_step[] = {durationPS_rosenbrock_CJ, durationPS_rkdopri};
		auto min_time_step = *min_element(time_step, time_step + 2);

		fprintf(file, "%s:\t\t<Total duration>\t\t<Time per step>\t\t<Number of steps> \n\n", name);
		fprintf(file, "\t\tRosenbrock (constant Jacobian):\t\t%.2g s (x%.2g)\t\t%.2g ms/step (x%.2g)\t\t%zu (x%.2g) \n",
				duration_rosenbrock_CJ, duration_rosenbrock_CJ / min_dur,
				durationPS_rosenbrock_CJ, durationPS_rosenbrock_CJ / min_time_step,
				n_steps_rosenbrock_CJ, (Double)n_steps_rosenbrock_CJ / min_steps);
		fprintf(file, "\t\tRKDOPRI:                       \t\t%.2g s (x%.2g)\t\t%.2g ms/step (x%.2g)\t\t%zu (x%.2g) \n\n",
				duration_rkdopri, duration_rkdopri / min_dur,
				durationPS_rkdopri, durationPS_rkdopri / min_time_step,
				n_steps_rkdopri, (Double)n_steps_rkdopri / min_steps);
		fclose(file);
	}
	else
	{
		/* Integration with Rosenbrock4 (using a constant Jacobian). */
		auto stepper_rosenbrock = make_dense_output(abs_err, rel_err, rosenbrock4<Double>());
		auto begin = steady_clock::now();
		for (size_t i = 0; i < n; ++i)
		{
			Y0[0] = i0;
			Y0[1] = a0;
			Y0[2] = m0;
			Y0[3] = s0;
			Y0[4] = z0;
			n_steps_rosenbrock_CJ = integrate_const(stepper_rosenbrock, make_pair(rhs(), jacobi_const()), Y0, Tstart, Tend, dt);
		}
		auto end = steady_clock::now();

		// Total time taken for n iterations in nanoseconds:
		auto duration_rosenbrock_CJ = duration_cast<duration<Double, nano>>(end - begin).count();
		// Time taken by one iteration in seconds:
		duration_rosenbrock_CJ = duration_rosenbrock_CJ / (1e9 * n);
		// Time taken by one step in milliseconds:
		Double durationPS_rosenbrock_CJ = (duration_rosenbrock_CJ / (Double)n_steps_rosenbrock_CJ) * 1e3;

		/* Integration with Rosenbrock4 (using a variable Jacobian). */
		begin = steady_clock::now();
		for (size_t i = 0; i < n; ++i)
		{
			Y0[0] = i0;
			Y0[1] = a0;
			Y0[2] = m0;
			Y0[3] = s0;
			Y0[4] = z0;
			n_steps_rosenbrock_VJ = integrate_const(stepper_rosenbrock, make_pair(rhs(), jacobi_var()), Y0, Tstart, Tend, dt);
		}
		end = steady_clock::now();

		// Total time taken for n iterations in nanoseconds:
		auto duration_rosenbrock_VJ = duration_cast<duration<Double, nano>>(end - begin).count();
		// Time taken by one iteration in seconds:
		duration_rosenbrock_VJ = duration_rosenbrock_VJ / (1e9 * n);
		// Time taken by one step in milliseconds:
		Double durationPS_rosenbrock_VJ = (duration_rosenbrock_VJ / (Double)n_steps_rosenbrock_VJ) * 1e3;

		/* Integration with Runge-Kutta-Dormand–Prince. */
		abs_err = 1.0e-12; // Reset absolute error.
		rel_err = 1.0e-12; // Reset relative error.

		auto stepper_rkdopri = make_dense_output(abs_err, rel_err, runge_kutta_dopri5<vector_type>());
		begin = steady_clock::now();
		for (size_t i = 0; i < n; ++i)
		{
			Y0[0] = i0;
			Y0[1] = a0;
			Y0[2] = m0;
			Y0[3] = s0;
			Y0[4] = z0;
			n_steps_rkdopri = integrate_const(stepper_rkdopri, rhs(), Y0, Tstart, Tend, dt);
		}
		end = steady_clock::now();

		// Total time taken for n iterations in nanoseconds:
		auto duration_rkdopri = duration_cast<duration<Double, nano>>(end - begin).count();
		// Time taken by one iteration in seconds:
		duration_rkdopri = duration_rkdopri / (1e9 * n);
		// Time taken by one step in milliseconds:
		Double durationPS_rkdopri = (duration_rkdopri / (Double)n_steps_rkdopri) * 1e3;

		/* Save results in file. */

		Double duration[] = {duration_rosenbrock_CJ, duration_rosenbrock_VJ, duration_rkdopri};
		auto min_dur = *min_element(duration, duration + 3);

		size_t steps[] = {n_steps_rosenbrock_CJ, n_steps_rosenbrock_VJ, n_steps_rkdopri};
		auto min_steps = (Double)*min_element(steps, steps + 3);

		Double time_step[] = {durationPS_rosenbrock_CJ, durationPS_rosenbrock_VJ, durationPS_rkdopri};
		auto min_time_step = *min_element(time_step, time_step + 3);

		fprintf(file, "%s:\t\t<Total duration>\t\t<Time per step>\t\t<Number of steps> \n\n", name);
		fprintf(file, "\tRosenbrock (constant Jacobian):\t\t%.2g s (x%.2g)\t\t%.2g ms/step (x%.2g)\t\t%zu (x%.2g) \n",
				duration_rosenbrock_CJ, duration_rosenbrock_CJ / min_dur,
				durationPS_rosenbrock_CJ, durationPS_rosenbrock_CJ / min_time_step,
				n_steps_rosenbrock_CJ, (Double)n_steps_rosenbrock_CJ / min_steps);
		fprintf(file, "\tRosenbrock (variable Jacobian):\t\t%.2g s (x%.2g)\t\t%.2g ms/step (x%.2g)\t\t%zu (x%.2g) \n",
				duration_rosenbrock_VJ, duration_rosenbrock_VJ / min_dur,
				durationPS_rosenbrock_VJ, durationPS_rosenbrock_VJ / min_time_step,
				n_steps_rosenbrock_VJ, (Double)n_steps_rosenbrock_VJ / min_steps);
		fprintf(file, "\tRKDOPRI:                       \t\t%.2g s (x%.2g)\t\t%.2g ms/step (x%.2g)\t\t%zu (x%.2g) \n\n",
				duration_rkdopri, duration_rkdopri / min_dur,
				durationPS_rkdopri, durationPS_rkdopri / min_time_step,
				n_steps_rkdopri, (Double)n_steps_rkdopri / min_steps);
		fclose(file);
	}
}

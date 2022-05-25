#include <cstdio>
#include <cmath>
#include <boost/numeric/odeint.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

using namespace std;
using namespace boost::numeric::odeint;

/*******************************************************************************************
* Model 008
*
* Alternative names: 'Model C', 'Self consistent Ascasibar et al. model (constant τS)'
*******************************************************************************************/

const char NAME[] = "Model 008"; // Model name

/* File initialization function */

FILE *file_ini(const char *name)
{
    FILE *file;
    file = fopen(name, "a");
    fprintf(file, "t \t i(t) \t a(t) \t m(t) \t z(t) \t s(t) \n");
    return file;
}

/* Witness function */

struct file_observer
{
    FILE *fp;

    file_observer(FILE *file) : fp{file} {}

    template <class State>
    void operator()(const State &y, double t) const
    {
        fprintf(
            fp,
            "%.15g \t %.15g \t %.15g \t %.15g \t %.15g \t %.15g \n",
            t, y[0], y[1], y[2], y[3], y[4]
        );
    }
};

/* Constants */

const double C1 = 3.03e-6; // [Gyr]
const double C2 = 8.74e-5; // [Gyr]
const double tau_S = 2.6;  // [Gyr]
const double Zsun = 0.0134;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.09;
const double eta_ion = 955.29;
const double eta_diss = 380.93;
const double R = 0.18;

/* Initial conditions */

const double g0 = 1e-3 * 14.8285; // [Mₒ pc^(-3)]
const double i0 = g0 * 0.6;       // [Mₒ pc^(-3)]
const double a0 = g0 * 0.2;       // [Mₒ pc^(-3)]
const double m0 = g0 * 0.2;       // [Mₒ pc^(-3)]
const double z0 = g0 * 1e-4;      // [Mₒ pc^(-3)]
const double s0 = g0 * 0.0;       // [Mₒ pc^(-3)]
vector_type Y0{5};

/* Parameters */

/* 
*	There are no parameters in this model.
*/

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
struct ode_system
{
    template <class State>
    void operator()(const State &y, State &dydt, double t)
    {
        (void)(t);

        /* Auxiliary equations */
        double g = y[0] + y[1] + y[2];
        double Z = y[3] / g;
        double tau_R = C1 / y[0];
        double tau_C = C2 / ((y[1] + y[2]) * (Z + Zeff));
        double recombination = y[0] / tau_R;
        double cloud_formation = y[1] / tau_C;
        double psi = y[2] / tau_S;

        /* ODE system */
        dydt[0] = -recombination + (eta_ion + R) * psi;
        dydt[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
        dydt[2] = cloud_formation - (1 + eta_diss) * psi;
        dydt[3] = (Zsn * R - Z) * psi;
        dydt[4] = (1 - R) * psi;
    }
};

/* Jacobian */

struct jacobian
{
    template <class State, class Matrix>
    void operator()(const State &y, Matrix &J, const double &t, State &dfdt)
    {
        (void)(t);

        J(0, 0) = -660066.0066 * y[0];
        J(0, 1) = 0;
        J(0, 2) = 367.4884615;
        J(0, 3) = 0;
        J(0, 4) = 0;

        J(1, 0) = 660066.0066 * y[0] + (11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2);
        J(1, 1) = (11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2) - 11441.6476 * y[1] * (0.0000134 + y[3] / (y[1] + y[0] + y[2])) - 11441.6476 * (y[1] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]));
        J(1, 2) = -220.9076923 + (11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2) - 11441.6476 * y[1] * (0.0000134 + y[3] / (y[1] + y[0] + y[2]));
        J(1, 3) = (-11441.6476 * y[1] * (y[1] + y[2])) / (y[1] + y[0] + y[2]);
        J(1, 4) = 0;

        J(2, 0) = (-11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2);
        J(2, 1) = (-11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2) + 11441.6476 * y[1] * (0.0000134 + y[3] / (y[1] + y[0] + y[2])) + 11441.6476 * (y[1] + y[2]) * (0.0000134 + y[3] / (y[1] + y[0] + y[2]));
        J(2, 2) = -146.8961538 - (11441.6476 * y[1] * (y[1] + y[2]) * y[3]) / pow(y[1] + y[0] + y[2], 2) + 11441.6476 * y[1] * (0.0000134 + y[3] / (y[1] + y[0] + y[2]));
        J(2, 3) = (11441.6476 * y[1] * (y[1] + y[2])) / (y[1] + y[0] + y[2]);
        J(2, 4) = 0;

        J(3, 0) = (0.3846153846 * y[2] * y[3]) / pow(y[1] + y[0] + y[2], 2);
        J(3, 1) = (0.3846153846 * y[2] * y[3]) / pow(y[1] + y[0] + y[2], 2);
        J(3, 2) = (0.3846153846 * y[2] * y[3]) / pow(y[1] + y[0] + y[2], 2) + 0.3846153846 * (0.0162 - (1. * y[3]) / (y[1] + y[0] + y[2]));
        J(3, 3) = (-0.3846153846 * y[2]) / (y[1] + y[0] + y[2]);
        J(3, 4) = 0;

        J(4, 0) = 0;
        J(4, 1) = 0;
        J(4, 2) = 0.3153846154;
        J(4, 3) = 0;
        J(4, 4) = 0;

        dfdt[0] = 0;
        dfdt[1] = 0;
        dfdt[2] = 0;
        dfdt[3] = 0;
        dfdt[4] = 0;
    }
};

/*******************************************************************************************
* Solving the system
*******************************************************************************************/

/* Integration span */

const double t_start = 0.0; // [Gyr]
const double t_end = 1.0;   // [Gyr]

int main()
{
    /* Initial time step */
    const double dt = (t_end - t_start) * 1e-5;

    /* Absolute tolerance */
    double abs_tol = 1.0e-6;

    /* Relative tolerance */
    double rel_tol = 1.0e-8;

    /* Initialization of the initial conditions */
    Y0[0] = i0;
    Y0[1] = a0;
    Y0[2] = m0;
    Y0[3] = z0;
    Y0[4] = s0;

    /* Integration with Rosenbrock4 */
    auto stepper_rosenbrock = make_dense_output(abs_tol, rel_tol, rosenbrock4<double>());

    /* Save data in file */
    FILE *file;
    file = file_ini("plots/cpp_rosenbrock4.dat");
    integrate_const(
        stepper_rosenbrock,
        make_pair(ode_system(), jacobian()),
        Y0,
        t_start,
        t_end,
        dt,
        file_observer(file)
    );
    fclose(file);

    printf("Stellar density after %gGyr (rosenbrock4): %.3g Mₒpc^(-2)\n", t_end, Y0[4]);

    /* Reset of the initial conditions */
    Y0[0] = i0;
    Y0[1] = a0;
    Y0[2] = m0;
    Y0[3] = z0;
    Y0[4] = s0;

    /* Integration with Runge-Kutta-Dormand–Prince */
    auto stepper_rkdopri = make_dense_output(
        abs_tol,
        rel_tol,
        runge_kutta_dopri5<vector_type>());

    /* Save data in file */
    file = file_ini("plots/cpp_rkdopri.dat");
    integrate_const(
        stepper_rkdopri,
        ode_system(),
        Y0,
        t_start,
        t_end,
        dt,
        file_observer(file)
    );
    fclose(file);

    printf("Stellar density after %gGyr (rkdopri): %.3g Mₒpc^(-3)", t_end, Y0[4]);
}

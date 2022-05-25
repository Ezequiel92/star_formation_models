#include <cstdio>
#include <cmath>
#include <boost/numeric/odeint.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

using namespace std;
using namespace boost::numeric::odeint;

/*******************************************************************************************
* Model 016
*
* Alternative names: 'Model G'
*******************************************************************************************/

const char NAME[] = "Model 016"; // Model name

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
            t, y[0], y[1], y[2], y[3], y[4]);
    }
};

/* Constants */

const double alpha = 1.0;
const double K = alpha * 810e-3; //[Gyr Mₒ^(1/2) pc^(-3/2)]
const double P = 40.4;           // [Mₒ^(-1) pc^(3) cm^(-3)]
const double Rh = 1.9;           // [cm^3 Gyr^-1]
const double SigmaNuB = 8158.0;  // [cm^3 Gyr^-1]
const double Zsun = 0.0134;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.09;
const double eta_ion = 955.29;
const double eta_diss = 380.93;
const double R = 0.18;

/* Initial conditions */

const double i0 = 0.6;
const double a0 = 0.2;
const double m0 = 0.2;
const double z0 = 1e-4;
const double s0 = 0.0;
vector_type Y0{5};

/* Parameters */

const double n = 1e-3; // [cm^(-3)]
const double Z = z0;
const double tot0 = 1e-3 / P; // [Mₒ^(-1) pc^(3)]

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
struct ode_system
{
    template <class State>
    void operator()(const State &y, State &dydt, double t)
    {
        (void)(t);

        /* Auxiliary equations */
        double tau_R = 1 / (n * SigmaNuB);
        double tau_C = Zsun / (2 * n * Rh * (Z + Zeff));
        double recombination = y[0] / tau_R;
        double cloud_formation = y[1] / tau_C;
        double tau_S = K / sqrt((y[0] + y[1] + y[2]) * tot0);
        double psi = y[2] / tau_S;

        /* ODE system */
        dydt[0] = -recombination + eta_ion * psi;
        dydt[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
        dydt[2] = cloud_formation - (1 + eta_diss) * psi;
        dydt[3] = (Zsn * R - Z) * psi;
        dydt[4] = psi;
    }
};

/* Jacobian */

struct jacobian
{
    template <class State, class Matrix>
    void operator()(const State &y, Matrix &J, const double &t, State &dfdt)
    {
        (void)(t);

        J(0, 0) = -8158. * n + (589.6851852 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(0, 1) = (589.6851852 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(0, 2) = (589.6851852 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2])) + 1179.37037 * sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(0, 3) = 0;
        J(0, 4) = 0;

        J(1, 0) = 8158. * n - (354.5432099 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(1, 1) = -283.5820896 * n * (0.0000134 + Z) - (354.5432099 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(1, 2) = (-354.5432099 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2])) - 709.0864198 * sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(1, 3) = 0;
        J(1, 4) = 0;

        J(2, 0) = (-235.7592593 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(2, 1) = 283.5820896 * n * (0.0000134 + Z) - (235.7592593 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(2, 2) = (-235.7592593 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2])) - 471.5185185 * sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(2, 3) = 0;
        J(2, 4) = 0;

        J(3, 0) = (0.6172839506 * tot0 * (0.0162 - 1. * Z) * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(3, 1) = (0.6172839506 * tot0 * (0.0162 - 1. * Z) * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(3, 2) = (0.6172839506 * tot0 * (0.0162 - 1. * Z) * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2])) + 1.234567901 * (0.0162 - 1. * Z) * sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(3, 3) = 0;
        J(3, 4) = 0;

        J(4, 0) = (0.6172839506 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(4, 1) = (0.6172839506 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2]));
        J(4, 2) = (0.6172839506 * tot0 * y[2]) / sqrt(tot0 * (y[1] + y[0] + y[2])) + 1.234567901 * sqrt(tot0 * (y[1] + y[0] + y[2]));
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
        file_observer(file));
    fclose(file);

    printf("Stellar fraction after %gGyr (rosenbrock4): %.3g\n", t_end, Y0[4]);

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
        file_observer(file));
    fclose(file);

    printf("Stellar fraction after %gGyr (rkdopri): %.3g", t_end, Y0[4]);
}

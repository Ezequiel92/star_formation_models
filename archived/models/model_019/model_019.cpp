#include <cstdio>
#include <cmath>
#include <boost/numeric/odeint.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

using namespace std;
using namespace boost::numeric::odeint;

/*******************************************************************************************
* Model 019
*******************************************************************************************/

const char NAME[] = "Model 019"; // Model name

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

const double C1 = 0.809;    // [Gyr Mₒ^(1/2) pc^(-3/2)]
const double C2 = 3.024e-6; // [Gyr Mₒ pc^(-3)]
const double C3 = 8.745e-5; // [Gyr Mₒ pc^(-3)]
const double Zsun = 0.0134;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.09;
const double eta_ion = 955.29;
const double eta_diss = 380.93;
const double R = 0.18;

/* Initial conditions */

const double i0 = 0.9;
const double a0 = 0.05;
const double m0 = 0.05;
const double z0 = 1e-4;
const double s0 = 0.0;
vector_type Y0{5};

/* Parameters */

const double rho0 = 0.9 / 40.46;
const double g0 = i0 + a0 + m0;

/* System of equations */

/*
    Ionized gas fraction:       i(t) / ρ -> y[0]
    Atomic gas fraction:        a(t) / ρ -> y[1]
    Molecular gas fraction:     m(t) / ρ -> y[2]
    Metal fraction:             z(t) / ρ -> y[3]
    Stellar fraction:           s(t) / ρ -> y[4]	

    where ρ = i(t) + a(t) + m(t) + s(t)
     
    Each equation has Gyr^(-1) [years^(-9)] as its units on the LHS and RHS
*/
struct ode_system
{
    template <class State>
    void operator()(const State &y, State &dydt, double t)
    {
        (void)(t);

        /* Auxiliary equations */
        double g = y[0] + y[1] + y[2];
        double tau_S = (C1 * g0) / sqrt(g * rho0);
        double tau_R = C2 / (y[0] * rho0);
        double tau_C = C3 / (g * rho0 * (y[3] + Zeff));
        double recombination = y[0] / tau_R;
        double cloud_formation = y[1] / tau_C;
        double psi = y[2] / tau_S;

        /* ODE system */
        dydt[0] = -recombination + (eta_ion + R) * psi;
        dydt[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
        dydt[2] = cloud_formation - (1 + eta_diss) * psi;
        dydt[3] = (Zsn * R - y[3]) * psi;
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

        J(0, 0) = -6.61e5 * rho0 * y[0] + (590. * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2])));
        J(0, 1) = (590. * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2])));
        J(0, 2) = (590. * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) + (1180. * sqrt(rho0 * (y[1] + y[0] + y[2]))) / g0;
        J(0, 3) = 0;
        J(0, 4) = 0;

        J(1, 0) = 6.61e5 * rho0 * y[0] - (354.7 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) - 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]);
        J(1, 1) = (-354.7 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) - 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]) - 1.143e4 * rho0 * (y[1] + y[0] + y[2]) * (0.0000134 + y[3]);
        J(1, 2) = (-354.7 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) - (709.3 * sqrt(rho0 * (y[1] + y[0] + y[2]))) / g0 - 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]);
        J(1, 3) = -1.143e4 * rho0 * y[1] * (y[1] + y[0] + y[2]);
        J(1, 4) = 0;

        J(2, 0) = (-235.8 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) + 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]);
        J(2, 1) = (-235.8 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) + 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]) + 1.143e4 * rho0 * (y[1] + y[0] + y[2]) * (0.0000134 + y[3]);
        J(2, 2) = (-235.8 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) - (471.7 * sqrt(rho0 * (y[1] + y[0] + y[2]))) / g0 + 1.143e4 * rho0 * y[1] * (0.0000134 + y[3]);
        J(2, 3) = 1.143e4 * rho0 * y[1] * (y[1] + y[0] + y[2]);
        J(2, 4) = 0;

        J(3, 0) = (0.6175 * rho0 * y[2] * (0.0162 - 1. * y[3])) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2])));
        J(3, 1) = (0.6175 * rho0 * y[2] * (0.0162 - 1. * y[3])) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2])));
        J(3, 2) = (0.6175 * rho0 * y[2] * (0.0162 - 1. * y[3])) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) + (1.235 * sqrt(rho0 * (y[1] + y[0] + y[2])) * (0.0162 - 1. * y[3])) / g0;
        J(3, 3) = (-1.235 * y[2] * sqrt(rho0 * (y[1] + y[0] + y[2]))) / g0;
        J(3, 4) = 0;

        J(4, 0) = (0.5064 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2])));
        J(4, 1) = (0.5064 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2])));
        J(4, 2) = (0.5064 * rho0 * y[2]) / (g0 * sqrt(rho0 * (y[1] + y[0] + y[2]))) + (1.013 * sqrt(rho0 * (y[1] + y[0] + y[2]))) / g0;
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

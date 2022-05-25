#include <cstdio>
#include <cmath>
#include <boost/numeric/odeint.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

using namespace std;
using namespace boost::numeric::odeint;

/*******************************************************************************************
* Model 012
*
* Alternative names: 'Full model (variable τS without infall and outflow)'
*******************************************************************************************/

const char NAME[] = "Model 012"; // Model name

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

const double K = 1 / 400.0;       // [Mₒ^(-2/3) pc^(4/3) Gyr^(-1)]
const double Twarm = 10000.0;     // [K]
const double Tcold = 100.0;       // [K]
const double T1 = 50000.0;        // [K]
const double Sigma_ion = 8e-4;    // [Mₒ pc^(-2)]
const double Sigma_diss = 1.5e-4; // [Mₒ pc^(-2)]
const double C2 = 0.074;          // [Mₒ^2 pc^(-4) Gyr]
const double C4 = 798e-3;         // [Mₒ^2 pc^(-4) Gyr]
const double alpha = 7.0;
const double eta_i_lim = 1000.0;
const double eta_d_lim = 100.0;
const double Zsun = 0.006;
const double Zeff = 1e-3 * Zsun;
const double Zsn = 0.06;
const double R = 0.17;

/* Initial conditions */

const double g0 = 1e-3 * 14.8285; // [Mₒ pc^(-2)]
const double i0 = g0 * 0.6;       // [Mₒ pc^(-2)]
const double a0 = g0 * 0.2;       // [Mₒ pc^(-2)]
const double m0 = g0 * 0.2;       // [Mₒ pc^(-2)]
const double z0 = g0 * 1e-4;      // [Mₒ pc^(-2)]
const double s0 = g0 * 0.0;       // [Mₒ pc^(-2)]
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

    Each equation has Mₒ pc^(-2) Gyr^(-1) [solar_mass * parsec^(-2) * years^(-9)]
	as its units on the LHS and RHS
*/
struct ode_system
{
    template <class State>
    void operator()(const State &y, State &dydt, double t)
    {
        /* Auxiliary equations */
        double g = y[0] + y[1] + y[2];
        double tot = g + y[4];
        double star_elem = y[2] + alpha * y[1];
        double tau_S = pow(star_elem, 1 / 3.0) / (K * pow(g, 1 / 3.0) * pow(tot, 2 / 3.0));
        double psi = star_elem / tau_S;
        double eta_ion = eta_i_lim * (1 - exp(-y[1] / Sigma_ion));
        double eta_diss = eta_d_lim * (1 - exp(-y[2] / Sigma_diss));
        double Z = y[3] / g;
        double tau_R = C2 * (1 + T1 * psi / (Twarm * g)) / (g * tot);
        double tau_C = C4 * (1 + T1 * psi / (Tcold * g)) * Zsun / (g * tot * (Z + Zeff));
        double recombination = y[0] / tau_R;
        double cloud_formation = y[1] / tau_C;

        /* ODE system */
        dydt[0] = -recombination + (eta_ion + R) * psi;
        dydt[1] = -cloud_formation + recombination + (eta_diss - eta_ion - alpha * y[1] / star_elem) * psi;
        dydt[2] = cloud_formation - (eta_diss + y[2] / star_elem) * psi;
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

        J(0, 0) = (0.001666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0008333333333 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (13.51351351 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (13.51351351 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(0, 1) = (0.001666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0008333333333 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.01166666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) + (3125. * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / exp(1250. * y[1]) + (13.51351351 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.05833333333 * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666)))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(0, 2) = (0.001666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0008333333333 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.001666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) + (13.51351351 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.008333333333 * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666)))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(0, 3) = 0;
        J(0, 4) = (0.001666666667 * (0.17 + 1000. * (1. - 1. / exp(1250. * y[1]))) * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.1126126126 * y[0] * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));

        J(1, 0) = (0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (13.51351351 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (13.51351351 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) + (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667)) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(1, 1) = (0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + 0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1.25e6 / exp(1250. * y[1]) + (49. * y[1]) / pow(7. * y[1] + y[2], 2) - 7. / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) + (0.0008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.01166666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) - (13.51351351 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.05833333333 * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666)))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) + (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (5.833333333 * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666))) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (208.855472 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(1, 2) = (0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + 0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (666666.6667 / exp(6666.666666667 * y[2]) + (7. * y[1]) / pow(7. * y[1] + y[2], 2)) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) + (0.0008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.001666666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) - (13.51351351 * y[0] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.008333333333 * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666)))) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (13.51351351 * y[0] * (y[1] + y[0] + y[2] + y[4])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) + (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.8333333333 * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666))) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) - (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(1, 3) = (-208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4])) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(1, 4) = (0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (-1000. * (1. - 1. / exp(1250. * y[1])) + 100. * (1. - 1. / exp(6666.666666667 * y[2])) - (7. * y[1]) / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) - (0.1126126126 * y[0] * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (13.51351351 * y[0] * (y[1] + y[0] + y[2])) / (1. + (0.0125 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (174.0462267 * y[1] * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) - (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));

        J(2, 0) = (-0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) - (0.0008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (100. * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) - (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667)) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(2, 1) = (-0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0175 * y[2] * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(7. * y[1] + y[2], 1.3333333333333333) - (0.0008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (100. * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (0.01166666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) - (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) - (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (5.833333333 * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666))) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472 * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(2, 2) = (-0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) - 0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (666666.6667 / exp(6666.666666667 * y[2]) - (1. * y[2]) / pow(7. * y[1] + y[2], 2) + 1 / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) - (0.0008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * (100. * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) - (0.001666666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (7. * y[1] + y[2])) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333) - (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * y[3]) / ((y[1] + y[0] + y[2]) * (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666))) - (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (y[1] + y[0] + y[2] + y[4]) * ((0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666)) / (pow(y[1] + y[0] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333)) - (0.8333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.8333333333 * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / (pow(7. * y[1] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2], 0.6666666666666666))) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666)) + (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(2, 3) = (208.855472 * y[1] * (y[1] + y[0] + y[2] + y[4])) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));
        J(2, 4) = (-0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (100. * (1. - 1. / exp(6666.666666667 * y[2])) + y[2] / (7. * y[1] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) - (174.0462267 * y[1] * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / pow(1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666), 2) + (208.855472 * y[1] * (y[1] + y[0] + y[2]) * (6.e-6 + y[3] / (y[1] + y[0] + y[2]))) / (1. + (1.25 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666));

        J(3, 0) = (0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * y[3]) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0102 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (0.0102 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2], 0.6666666666666666);
        J(3, 1) = (0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * y[3]) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0102 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (0.0102 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.01166666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (0.0102 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(7. * y[1] + y[2], 0.3333333333333333);
        J(3, 2) = (0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * y[3]) / pow(y[1] + y[0] + y[2], 1.6666666666666667) + (0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0102 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0008333333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (0.0102 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.001666666667 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666) * (0.0102 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(7. * y[1] + y[2], 0.3333333333333333);
        J(3, 3) = (-0.0025 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666);
        J(3, 4) = (0.001666666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333) * (0.0102 - (1. * y[3]) / (y[1] + y[0] + y[2]))) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333);

        J(4, 0) = (0.001383333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0006916666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666);
        J(4, 1) = (0.001383333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0006916666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.009683333333 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333);
        J(4, 2) = (0.001383333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333) + (0.0006916666667 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(y[1] + y[0] + y[2], 0.6666666666666666) + (0.001383333333 * pow(y[1] + y[0] + y[2], 0.3333333333333333) * pow(y[1] + y[0] + y[2] + y[4], 0.6666666666666666)) / pow(7. * y[1] + y[2], 0.3333333333333333);
        J(4, 3) = 0;
        J(4, 4) = (0.001383333333 * pow(7. * y[1] + y[2], 0.6666666666666666) * pow(y[1] + y[0] + y[2], 0.3333333333333333)) / pow(y[1] + y[0] + y[2] + y[4], 0.3333333333333333);

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
        runge_kutta_dopri5<vector_type>()
    );

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

    printf("Stellar density after %gGyr (rkdopri): %.3g Mₒpc^(-2)", t_end, Y0[4]);
}

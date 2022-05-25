/*!
 * \file        src/star_formation/ez_sfr/ez_sfr.c
 * \date        03/2022
 * \brief       Compute the probability of star formation for a given gas cell.
 * \details     This file contains the routines to compute the stellar fraction that would form
 *              in a gas cell, according to our star formation model. It evolves a set of ODEs,
 *              that model the mass exchange between different phases of Hydrogen, metals, and stars.
 *              contains functions:
 *                double **read_ftable(const char *filename, size_t *rows, size_t *cols)
 *                void check_error(const int err_type, const int err_code)
 *                double eval_interp(struct InterpFunc *interp, const double x)
 *                void get_eta(glob_t globbuf, const double age, const int ion, struct InterpFunc *interp_func)
 *                int sf_ode(double t, const double y[], double f[], void *params)
 *                int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params)
 *                void fractions(const double int_time, struct Parameters *params, double y[])
 *                double stellar_fraction(const int index, const double int_time)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

// #include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <glob.h>
// #include <stdbool.h>
// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_spline.h>
// #include <gsl/gsl_interp.h>
// #include <gsl/gsl_matrix.h>
// #include <gsl/gsl_odeiv2.h>
// #include <gsl/gsl_integration.h>

// #include "../allvars.h"
#include "ez_sfr.h"

// #ifdef EZ_SFR

/*! \brief This function reads a file and saves its contents in a matrix.
 *
 *  Read a text file with space separated values, and stores the numerical
 *  data as doubles in a matrix. Each line in the file must be at most 1024
 *  characters long.
 *
 *  \param[in] filename Path to the file.
 *  \param[out] rows Where the number of rows will be stored.
 *  \param[out] cols Where the number of columns will be stored.
 *
 *  \return A pointer to the matrix storing the data, or NULL if the request failed.
 */
double **read_ftable(const char *filename, size_t *rows, size_t *cols)
{
    if (rows == NULL || cols == NULL || filename == NULL)
        return NULL;

    *rows = 0;
    *cols = 0;

    FILE *fp = fopen(filename, "r");

    if (fp == NULL)
    {
        fprintf(stdout, "Could not open %s\n", filename);
        return NULL;
    }

    double **matrix = NULL, **tmp;
    char line[1024];

    while (fgets(line, sizeof(line), fp))
    {
        if (*cols == 0)
        {
            /* Determine the number of columns based on the first row */
            char *scan = line;
            double dummy;
            int offset = 0;
            while (sscanf(scan, "%lf%n", &dummy, &offset) == 1)
            {
                scan += offset;
                (*cols)++;
            }
        }

        tmp = (double **)realloc(matrix, (*rows + 1) * sizeof(*matrix));

        if (tmp == NULL)
        {
            fclose(fp);
            fprintf(stdout, "The allocation of memory for a matrix failed\n");
            return matrix;
        }

        matrix = tmp;

        matrix[*rows] = (double *)calloc(*cols, sizeof(*matrix[*rows]));

        if (matrix[*rows] == NULL)
        {
            fclose(fp);
            fprintf(stdout, "The allocation of memory for a matrix failed\n");
            if (*rows == 0)
            {
                free(matrix);
                return NULL;
            }

            return matrix;
        }

        int offset = 0;
        char *scan = line;

        for (size_t i = 0; i < *cols; ++i)
        {
            if (sscanf(scan, "%lf%n", matrix[*rows] + i, &offset) == 1)
                scan += offset;
            else
                matrix[*rows][i] = 0; /* Set element that could not be read to 0 */
        }

        (*rows)++;
    }

    fclose(fp);

    return matrix;
}

/*! \brief Prints a human redable explanation of known error codes.
 *
 *  Knows about errors from GSL and GLOB.
 *
 *  \param[in] err_type Code that indicates the library from which the error originates.
 *  \param[in] err_code error code.
 *
 *  \return void.
 */
void check_error(const int err_type, const int err_code)
{
    switch (err_type)
    {
    case GLOB_ERROR:
        switch (err_code)
        {
        case GLOB_NOSPACE:
            printf("GLOB_NOSPACE\n");
            printf("Glob run out of memory\n");
            break;
        case GLOB_ABORTED:
            printf("GLOB_ABORTED\n");
            printf("Glob had a read error\n");
            break;
        case GLOB_NOMATCH:
            printf("GLOB_NOMATCH\n");
            printf("Glob found no matches\n");
            break;
        default:
            printf("Unknown GLOB error\n");
            break;
        }
        break;
    case GLS_ODE_ERROR:
        switch (err_code)
        {
        case GSL_EBADFUNC:
            printf("GSL_EBADFUNC\n");
            printf("`sf_ode()` or `jacobian()` could not be executed successfully\n");
            break;
        case GSL_FAILURE:
            printf("GSL_FAILURE\n");
            printf("The stepping function was unable to compute the requested step (required step may be too small)\n");
            break;
        case GSL_EFAULT:
            printf("GSL_EFAULT\n");
            printf("Driver object is not appropriately set\n");
            break;
        case GSL_EMAXITER:
            printf("GSL_EMAXITER\n");
            printf("Maximum number of steps has been reached\n");
            break;
        case GSL_ENOPROG:
            printf("GSL_ENOPROG\n");
            printf("The step size has dropped below its minimum value\n");
            break;
        default:
            printf("Unknown GSL ODE solver error\n");
            break;
        }
        break;
    case GSL_INT_ERROR:
        switch (err_code)
        {
        case GSL_EMAXITER:
            printf("GSL_EMAXITER\n");
            printf("The maximum number of subdivisions was exceeded\n");
            break;
        case GSL_EROUND:
            printf("GSL_EROUND\n");
            printf("Cannot reach tolerance because of roundoff error\n");
            break;
        case GSL_ESING:
            printf("GSL_ESING\n");
            printf("A non-integrable singularity or other bad integrand behavior was found in the integration interval\n");
            break;
        case GSL_EDIVERGE:
            printf("GSL_EDIVERGE\n");
            printf("The integral is divergent, or too slowly convergent to be integrated numerically\n");
            break;
        case GSL_EDOM:
            printf("GSL_EDOM\n");
            printf("Error in the values of the input arguments\n");
            break;
        default:
            printf("Unknown GSL integration error\n");
            break;
        }
        break;
    default:
        printf("Unknown error\n");
        break;
    }
    printf("Error return value = %d\n", err_code);
    fflush(stdout);
    endrun(8888);
}

/*! \brief Evaluate a given interpolation function.
 *
 *  Evaluate a interpolation function at a given value.
 *  If the value is out of range, the function at its
 *  extrema is used.
 *
 *  \param[in] interp Interpolation function.
 *  \param[in] x Value at which the interpolation function will be evaluated.
 *
 *  \return The result of evaluating the interpolation function.
 */
double eval_interp(struct InterpFunc *interp, const double x)
{
    if (isnan(interp->x_min) || isnan(interp->x_max))
    {
        fprintf(stdout, "x_min and/or x_max where not stored correctly\n");
        fflush(stdout);
        endrun(8888);
    }
    else if (x < interp->x_min)
    {
        return gsl_spline_eval(interp->interp, interp->x_min, interp->acc);
    }
    else if (x > interp->x_max)
    {
        return gsl_spline_eval(interp->interp, interp->x_max, interp->acc);
    }
    else
    {
        return gsl_spline_eval(interp->interp, x, interp->acc);
    }
}

/*! \brief Compute an interpolation function for a photodissociation efficiency.
 *
 *  Construct the interpolation function for the parameter `eta` as a function of the metallicity.
 *  There are two possible `eta`s, `eta_ion`: ionization efficiency for Hydrogen atoms,
 *  and `eta_diss`: photodissociation efficiency for Hydrogen molecules.
 *
 *  \param[in] globbuf Buffer with the paths to the files in the folder ./Q_data.
 *  \param[in] age Maximum age for the stellar population, in Gyr.
 *  \param[in] ion 0: Compute the ionization efficiency for Hydrogen atoms (`eta_ion`),
                   1: Compute the photodissociation efficiency for Hydrogen molecules (`eta_diss`).
 *  \param[out] interp_func Where the resulting interpolation function will be stored.
 *
 *  \return void.
 */
void get_eta(glob_t globbuf, const double age, const int ion, struct InterpFunc *interp_func)
{
    size_t cols, rows, size, j = 0;
    double **eta_table, etas[globbuf.gl_pathc], metallicities[globbuf.gl_pathc], log_age = log10(age * 1.0e9);
    char **paths = globbuf.gl_pathv;
    struct InterpFunc eta_age_interp;

    while (*paths)
    {
        eta_table = read_ftable(*paths, &rows, &cols);

        if (eta_table == NULL)
        {
            fprintf(stdout, "Could not read the files with the Q values\n");
            fflush(stdout);
            endrun(8888);
        }

        metallicities[j] = atof(basename(*paths));

        /* Calculate the number of rows with ages below `age` */
        size = 0;
        while (size < rows && eta_table[size][0] <= log_age)
            size++;

        double age_col[size];
        double eta_col[size];

        for (size_t i = 0; i < size; ++i)
        {
            age_col[i] = eta_table[i][0];
            eta_col[i] = eta_table[i][ion + 1];
        }

        eta_age_interp.acc = gsl_interp_accel_alloc();
        eta_age_interp.interp = gsl_spline_alloc(gsl_interp_linear, size);
        eta_age_interp.x_min = eta_table[0][0];
        eta_age_interp.x_max = eta_table[rows - 1][0];

        gsl_spline_init(eta_age_interp.interp, age_col, eta_col, size);

        etas[j] = eval_interp(&eta_age_interp, log_age);

        for (size_t i = 0; i < rows; ++i)
            free(eta_table[i]);
        free(eta_table);
        gsl_interp_accel_free(eta_age_interp.acc);
        gsl_spline_free(eta_age_interp.interp);

        paths++;
        j++;
    }

    gsl_spline_init(interp_func->interp, metallicities, etas, globbuf.gl_pathc);
    interp_func->x_min = metallicities[0];
    interp_func->x_max = metallicities[j - 1];
}

/*! \brief Evaluate the systems of equations of our star formation model.
 *
 *  Evaluate the five ODEs of our model, using the following variables:
 *
 *  Ionized gas fraction:       i(t) / ρ -> y[0]
 *  Atomic gas fraction:        a(t) / ρ -> y[1]
 *  Molecular gas fraction:     m(t) / ρ -> y[2]
 *  Metal fraction:             z(t) / ρ -> y[3]
 *  Stellar fraction:           s(t) / ρ -> y[4]
 *
 *  where ρ = i(t) + a(t) + m(t) + s(t), and each equation has units of Gyr^(-1).
 *
 *  The parameters are `rho0` [Mₒ pc^(-3)], `g0` [dimensionless], and the interpolation
 *  functions for `eta_ion` and `eta_diss`.
 *
 *  \param[in] t Unused variable to conform to the gsl_odeiv2_driver_alloc_y_new() API.
 *  \param[in] y Values of the variables at which the ODEs will be evaluated.
 *  \param[out] f Where the results of evaluating the ODEs will be stored.
 *  \param[in] params Constant parameters for the ODEs.
 *
 *  \return Constant GSL_SUCCESS to confirm that the computation was successful.
 */
int sf_ode(double t, const double y[], double f[], void *params)
{
    (void)(t);

    struct Parameters *parameters = (struct Parameters *)params;

    double rho0 = parameters->rho0; /* Mean total density */
    double g0 = parameters->g0;     /* Initial gas fraction: (i(0) + a(0) + m(0)) / ρ */
    struct InterpFunc *interp_ion = parameters->interp_ion;
    struct InterpFunc *interp_diss = parameters->interp_diss;

    /* Auxiliary equations */
    double eta_ion = eval_interp(interp_ion, y[3]);
    double eta_diss = eval_interp(interp_diss, y[3]);
    double g = y[0] + y[1] + y[2];
    double tau_S = (C1 * g0) / sqrt(g * rho0);
    double tau_R = C2 / (y[0] * rho0);
    double tau_C = C3 / (g * rho0 * (y[3] + Zeff));
    double recombination = y[0] / tau_R;
    double cloud_formation = y[1] / tau_C;
    double psi = y[2] / tau_S;

    /* ODE system */
    f[0] = -recombination + (eta_ion + R) * psi;
    f[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
    f[2] = cloud_formation - (1 + eta_diss) * psi;
    f[3] = (Zsn * R - y[3]) * psi;
    f[4] = (1 - R) * psi;

    return GSL_SUCCESS;
};

/*! \brief Evaluate the jacobian for our star formation model.
 *
 *  Evaluate the jacobian matrix of our model, using the following variables:
 *
 *  Ionized gas fraction:       i(t) / ρ -> y[0]
 *  Atomic gas fraction:        a(t) / ρ -> y[1]
 *  Molecular gas fraction:     m(t) / ρ -> y[2]
 *  Metal fraction:             z(t) / ρ -> y[3]
 *  Stellar fraction:           s(t) / ρ -> y[4]
 *
 *  where ρ = i(t) + a(t) + m(t) + s(t).
 *
 *  The parameters are `rho0` [Mₒ pc^(-3)], `g0` [dimensionless], and the interpolation
 *  functions for `eta_ion` and `eta_diss`.
 *
 *  \param[in] t Unused variable to conform to the gsl_odeiv2_driver_alloc_y_new() API.
 *  \param[in] y Values of the variables at which the jacobian will be evaluated.
 *  \param[out] dfdy Where the results of evaluating the jacobian will be stored.
 *  \param[out] dfdt Where the result of evaluating the time derivatives will be stored.
 *  \param[in] params Constant parameters for the jacobian.
 *
 *  \return Constant GSL_SUCCESS to confirm that the computation was successful.
 */
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    (void)(t);

    struct Parameters *parameters = (struct Parameters *)params;

    double rho0 = parameters->rho0;
    double g0 = parameters->g0;
    struct InterpFunc *interp_ion = parameters->interp_ion;
    struct InterpFunc *interp_diss = parameters->interp_diss;

    double eta_ion = eval_interp(interp_ion, y[3]);
    double eta_diss = eval_interp(interp_diss, y[3]);

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 5, 5);
    gsl_matrix *m = &dfdy_mat.matrix;

    gsl_matrix_set(m, 0, 0, (0.049383732822699825 * rho0 * y[2] + 0.5 * rho0 * eta_ion * y[2] - 535200.887413877 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 0, 1, ((0.049383732822699825 + 0.5 * eta_ion) * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 0, 2, (0.09876746564539965 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) + pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) * eta_ion + 0.049383732822699825 * rho0 * y[2] + 0.5 * rho0 * eta_ion * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 0, 3, 0);
    gsl_matrix_set(m, 0, 4, 0);

    gsl_matrix_set(m, 1, 0, (-0.5 * rho0 * eta_ion * y[2] + 0.5 * rho0 * eta_diss * y[2] + 535200.887413877 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0] - 0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] - 9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 1, 1, (-0.5 * rho0 * eta_ion * y[2] + 0.5 * rho0 * eta_diss * y[2] - 0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0] - 0.2479695231261206 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] - 0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[2] - 9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0] * y[3] - 18505.188292994077 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3] - 9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[2] * y[3]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 1, 2, (pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) * eta_diss - 1.0 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) * eta_ion - 0.5 * rho0 * eta_ion * y[2] + 0.5 * rho0 * eta_diss * y[2] - 0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] - 9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 1, 3, (-11435.019367644041 * y[0] - 11435.019367644041 * y[1] - 11435.019367644041 * y[2]) * rho0 * y[1]);
    gsl_matrix_set(m, 1, 4, 0);

    gsl_matrix_set(m, 2, 0, (-0.5 * rho0 * y[2] - 0.5 * rho0 * eta_diss * y[2] + 0.12398476156306032 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] + 9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 2, 1, (-0.5 * rho0 * y[2] - 0.5 * rho0 * eta_diss * y[2] + 0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0] + 0.2479695231261206 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] + 0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[2] + 9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0] * y[3] + 18505.188292994077 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3] + 9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[2] * y[3]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 2, 2, (-1.0 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) - 1.0 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) * eta_diss - 0.5 * rho0 * y[2] - 0.5 * rho0 * eta_diss * y[2] + 0.12398476156306032 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] + 9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 2, 3, (11435.019367644041 * y[0] + 11435.019367644041 * y[1] + 11435.019367644041 * y[2]) * rho0 * y[1]);
    gsl_matrix_set(m, 2, 4, 0);

    gsl_matrix_set(m, 3, 0, ((0.007156218517396961 - 0.5 * y[3]) * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 3, 1, ((0.007156218517396961 - 0.5 * y[3]) * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 3, 2, (0.014312437034793922 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) - 1.0 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) * y[3] + 0.007156218517396961 * rho0 * y[2] - 0.5 * rho0 * y[2] * y[3]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 3, 3, (-1 * sqrt((y[0] + y[1] + y[2]) * rho0) * y[2]) / (0.8091454722567165 * g0));
    gsl_matrix_set(m, 3, 4, 0);

    gsl_matrix_set(m, 4, 0, (0.4506162671773002 * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 4, 1, (0.4506162671773002 * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 4, 2, (0.9012325343546004 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) + 0.4506162671773002 * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
    gsl_matrix_set(m, 4, 3, 0);
    gsl_matrix_set(m, 4, 4, 0);

    dfdt[0] = 0;
    dfdt[1] = 0;
    dfdt[2] = 0;
    dfdt[3] = 0;
    dfdt[4] = 0;

    return GSL_SUCCESS;
};

/*! \brief Numerical integration of the ODEs.
 *
 *  Integrate the system of ODEs, to obtain the following fractions:
 *
 *  Ionized gas fraction:       i(t) / ρ -> y[0]
 *  Atomic gas fraction:        a(t) / ρ -> y[1]
 *  Molecular gas fraction:     m(t) / ρ -> y[2]
 *  Metal fraction:             z(t) / ρ -> y[3]
 *  Stellar fraction:           s(t) / ρ -> y[4]
 *
 *  where ρ = i(t) + a(t) + m(t) + s(t).
 *
 *  \param[in] int_time Integration time, in Gyr.
 *  \param[in] params ODE parameters.
 *  \param[out] y Where the final fraction of each phase will be stored.
 *
 *  \return void.
 */
void fractions(const double int_time, struct Parameters *params, double y[])
{
    gsl_odeiv2_system sys = {sf_ode, jacobian, NUMEQU, params};
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf, int_time * 1e-5, abs_tol, rel_tol);

    /* Solve the ODE system using the BDF method */
    double t0 = 0.0;
    int gsl_err = gsl_odeiv2_driver_apply(driver, &t0, int_time, y);

    if (gsl_err != GSL_SUCCESS)
    {
        printf("GSL ODE solver error: ");
        check_error(GLS_ODE_ERROR, gsl_err);
    }

    gsl_odeiv2_driver_free(driver);
};

/*! \brief Compute the stellar fraction according to our model.
 *
 *  This fraction can be used to calculate the probability of forming a star population.
 *
 *  \param[in] index Index of the gas cell in question.
 *  \param[in] int_time Integration time in internal units.
 *
 *  \return The final stellar fraction ∈ [0, 1].
 */
double stellar_fraction(const int index, const double int_time)
{
    /* T [internal_units] * t_s = T [s] */
    const double t_s = All.UnitTime_in_s / All.HubbleParam;
    /* T [internal_units] * t_Gyr = T [Gyr] */
    const double t_Gyr = t_s / (SEC_PER_YEAR * 1e9);
    /* RHO [internal_units] * rho_cgs = RHO [g cm^(-3)] */
    const double rho_cgs = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;
    /* RHO [internal_units] * rho_cosmo = RHO [Mₒ pc^(-3)] */
    const double rho_cosmo = rho_cgs * PARSEC * PARSEC * PARSEC / SOLAR_MASS;

    double y[NUMEQU];
    double nh0 = 0, coolrate = 0;
    double ne = SphP[index].Ne;
    // SetOutputGasState(index, &ne, &nh0, &coolrate);

    /* Total baryonic density (tot0 = ρ_i(0) + ρ_a(0) + ρ_m(0) + ρ_s(0)) in [Mₒ pc^(-3)] */
    double tot0 = SphP[index].d.Density * rho_cosmo;

    double star_fraction = 0.0;

    glob_t globbuf;
    int glob_err = glob(IMF_PATH, GLOB_ERR, NULL, &globbuf);

    if (glob_err != 0)
    {
        printf("GLOB error: ");
        check_error(GLOB_ERROR, glob_err);
    }

    struct InterpFunc interp_ion;
    struct InterpFunc interp_diss;

    interp_ion.acc = gsl_interp_accel_alloc();
    interp_ion.interp = gsl_spline_alloc(gsl_interp_linear, globbuf.gl_pathc);
    interp_ion.x_min = NAN;
    interp_ion.x_max = NAN;

    get_eta(globbuf, int_time * t_Gyr, 0, &interp_ion);

    interp_diss.acc = gsl_interp_accel_alloc();
    interp_diss.interp = gsl_spline_alloc(gsl_interp_linear, globbuf.gl_pathc);
    interp_diss.x_min = NAN;
    interp_diss.x_max = NAN;
    get_eta(globbuf, int_time * t_Gyr, 1, &interp_diss);

    struct Parameters params;

    params.g0 = 1.0;
    params.interp_ion = &interp_ion;
    params.interp_diss = &interp_diss;

    for (int i = 0; i < DIVISIONS; ++i)
    {
        y[0] = SphP[0].fHII; /* Ionized fraction */
        y[1] = SphP[0].fHI;  /* Atomic fraction */
        y[2] = SphP[0].fmol; /* Molecular fraction */
        y[3] = SphP[0].fZ;   /* Metal fraction */
        y[4] = SphP[0].fstar;
        params.rho0 = tot0 * f_rho[i];
        fractions(int_time * t_Gyr, &params, y);
        star_fraction += y[4] * rho_pdf[i];
    }

    gsl_interp_accel_free(interp_ion.acc);
    gsl_spline_free(interp_ion.interp);
    gsl_interp_accel_free(interp_diss.acc);
    gsl_spline_free(interp_diss.interp);

    return star_fraction; /* Return final stellar mass fraction, formed after `int_time` [Gyr] */
}

// #endif /* of EZ_SFR */
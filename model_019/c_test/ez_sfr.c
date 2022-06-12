/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/ez_sfr/ez_sfr.c
 * \date        06/2022
 * \brief       Compute the star formation rate for a given gas cell.
 * \details     This file contains the routines to compute the star formation rate, according to our
 *              star formation model. It evolves a set of ODEs, that model the mass exchange between
 *              different phases of Hydrogen, metals, and stars.
 *              contains functions:
 *                double **read_ftable(const char *filename, size_t *rows, size_t *cols)
 *                void check_error(const int err_type, const int err_code)
 *                double eval_interp(struct InterpFunc *interp, const double x)
 *                void get_eta(glob_t globbuf, const double age, const int ion, struct InterpFunc *interp_func)
 *                int sf_ode(double t, const double y[], double f[], void *ode_params)
 *                int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *ode_params)
 *                void fractions(const double int_time, struct ODEParameters *ode_params, double y[])
 *                double rate_of_star_formation(const int index)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "./ez_sfr.h"

// #ifdef EZ_SFR

/*! \brief This function loads a text file into memory as a C matrix.
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
  if(rows == NULL || cols == NULL || filename == NULL)
    {
      return NULL;
    }

  *rows = 0;
  *cols = 0;

  FILE *fp = fopen(filename, "r");

  if(fp == NULL)
    {
      printf("Could not open %s\n", filename);
      return NULL;
    }

  double **matrix = NULL, **tmp;
  char line[1024];

  while(fgets(line, sizeof(line), fp))
    {
      if(*cols == 0)
        {
          /* Determine the number of columns based on the first row */
          char *scan = line;
          double dummy;
          int offset = 0;
          while(sscanf(scan, "%lf%n", &dummy, &offset) == 1)
            {
              scan += offset;
              (*cols)++;
            }
        }

      tmp = (double **)realloc(matrix, (*rows + 1) * sizeof(*matrix));

      if(tmp == NULL)
        {
          fclose(fp);
          printf("The allocation of memory in read_ftable() failed\n");
          return matrix;
        }

      matrix = tmp;

      matrix[*rows] = (double *)calloc(*cols, sizeof(*matrix[*rows]));

      if(matrix[*rows] == NULL)
        {
          fclose(fp);
          printf("The allocation of memory in read_ftable() failed\n");
          if(*rows == 0)
            {
              free(matrix);
              return NULL;
            }

          return matrix;
        }

      int offset = 0;
      char *scan = line;

      for(size_t i = 0; i < *cols; ++i)
        {
          if(sscanf(scan, "%lf%n", matrix[*rows] + i, &offset) == 1)
            {
              scan += offset;
            }
          else
            {
              matrix[*rows][i] = 0; /* Set elements that could not be read to 0 */
            }
        }

      (*rows)++;
    }

  fclose(fp);

  return matrix;
}

/*! \brief Prints a human readable explanation of some error codes.
 *
 *  Knows about errors in GSL and GLOB.
 *
 *  \param[in] err_type Number that indicates the library from which the error originates.
 *  \param[in] err_code Error code.
 *
 *  \return void.
 */
void check_error(const int err_type, const int err_code)
{
  switch(err_type)
    {
      case GLOB_ERROR:
        switch(err_code)
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
        switch(err_code)
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
      case GSL_INTERP_ERROR:
        switch(err_code)
          {
            case GSL_EDOM:
              printf("GSL_EDOM\n");
              printf("The input argument is outside the interpolation range\n");
              break;
            default:
              printf("Unknown GSL interpolation error\n");
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
 *  If the value is out of range, the function at its extrema is used.
 *
 *  \param[in] interp Interpolation function.
 *  \param[in] x Value at which the interpolation function will be evaluated.
 *
 *  \return The result of evaluating the interpolation function.
 */
double eval_interp(struct InterpFunc *interp, const double x)
{
  double y = 0;

  if(isnan(interp->x_min) || isnan(interp->x_max))
    {
      fprintf(stdout, "x_min and/or x_max in the InterpFunc struct where not stored correctly\n");
      fflush(stdout);
      endrun(8888);
    }
  else if(x < interp->x_min)
    {
      y = gsl_spline_eval(interp->interp, interp->x_min, interp->acc);
    }
  else if(x > interp->x_max)
    {
      y = gsl_spline_eval(interp->interp, interp->x_max, interp->acc);
    }
  else
    {
      y = gsl_spline_eval(interp->interp, x, interp->acc);
    }

  //if(y == GSL_NAN)
    //{
      //printf("GSL interpolation error: ");
      //check_error(GSL_INTERP_ERROR, y);
    //}

  return y;
}

/*! \brief Compute an interpolation function for the photodissociation efficiencies.
 *
 *  Compute the interpolation function for the parameter `eta` as a function of the metallicity.
 *  There are two possible `eta`s:
 *  `eta_ion`: ionization efficiency for Hydrogen atoms.
 *  `eta_diss`: photodissociation efficiency for Hydrogen molecules.
 *
 *  \param[in] globbuf Buffer with the paths to the files in the folder ./eta_data.
 *  \param[in] age Maximum age for the stellar population, in Gyr.
 *  \param[in] ion 0: Compute the ionization efficiency for Hydrogen atoms (`eta_ion`),
                   1: Compute the photodissociation efficiency for Hydrogen molecules (`eta_diss`).
 *  \param[out] interp_func Where the resulting interpolation function will be stored.
 *
 *  \return void.
 */
void get_eta(glob_t globbuf, const double age, const int ion, struct InterpFunc *interp_func)
{
  size_t cols, rows, j;
  double **eta_table, etas[globbuf.gl_pathc], metallicities[globbuf.gl_pathc], log_age = log10(age * 1.0e9);
  char **paths = globbuf.gl_pathv;
  struct InterpFunc eta_age_interp;

  j = 0;
  while(*paths)
    {
      eta_table = read_ftable(*paths, &rows, &cols);

      if(eta_table == NULL || cols != 3)
        {
          fprintf(stdout, "Something went wrong while trying to read the eta files\n");
          fflush(stdout);
          endrun(8888);
        }

      metallicities[j] = atof(basename(*paths));

      double age_col[rows];
      double eta_col[rows];

      for(size_t i = 0; i < rows; ++i)
        {
          age_col[i] = eta_table[i][0];
          eta_col[i] = eta_table[i][ion + 1];
        }

      eta_age_interp.acc    = gsl_interp_accel_alloc();
      eta_age_interp.interp = gsl_spline_alloc(gsl_interp_linear, rows);
      eta_age_interp.x_min  = eta_table[0][0];
      eta_age_interp.x_max  = eta_table[rows - 1][0];

      gsl_spline_init(eta_age_interp.interp, age_col, eta_col, rows);

      etas[j] = eval_interp(&eta_age_interp, log_age);

      for(size_t i = 0; i < rows; ++i)
        {
          free(eta_table[i]);
        }
      free(eta_table);
      gsl_spline_free(eta_age_interp.interp);
      gsl_interp_accel_free(eta_age_interp.acc);

      paths++;
      j++;
    }

  gsl_spline_init(interp_func->interp, metallicities, etas, globbuf.gl_pathc);
  interp_func->x_min = metallicities[0];
  interp_func->x_max = metallicities[globbuf.gl_pathc - 1];
}

/*! \brief Evaluate the systems of equations for the star formation model.
 *
 *  Evaluate the five ODEs of the model, using the following variables:
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
 *  \param[in] ode_params Constant parameters for the ODEs.
 *
 *  \return Constant GSL_SUCCESS to confirm that the computation was successful.
 */
int sf_ode(double t, const double y[], double f[], void *ode_params)
{
  (void)(t);

  struct ODEParameters parameters = *(struct ODEParameters *)ode_params;
  double rho0                     = parameters.rho0; /* Mean total density */
  double g0                       = parameters.g0;   /* Initial gas fraction: (i(0) + a(0) + m(0)) / ρ */
  double eta_ion                  = eval_interp(parameters.interp_ion, y[3]);
  double eta_diss                 = eval_interp(parameters.interp_diss, y[3]);

  /* Auxiliary equations */
  double g               = y[0] + y[1] + y[2];
  double tau_S           = (const_1 * g0) / sqrt(g * rho0);
  double tau_R           = const_2 / (y[0] * rho0);
  double tau_C           = const_3 / (g * rho0 * (y[3] + Zeff));
  double recombination   = y[0] / tau_R;
  double cloud_formation = y[1] / tau_C;
  double psi             = y[2] / tau_S;

  /* ODE system */
  f[0] = -recombination + (eta_ion + Rf) * psi;
  f[1] = -cloud_formation + recombination + (eta_diss - eta_ion) * psi;
  f[2] = cloud_formation - (1 + eta_diss) * psi;
  f[3] = (Zsn * Rf - y[3]) * psi;
  f[4] = (1 - Rf) * psi;

  return GSL_SUCCESS;
};

/*! \brief Evaluate the jacobian of the star formation model.
 *
 *  Evaluate the jacobian matrix of the model, using the following variables:
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
 *  \param[in] ode_params Constant parameters for the jacobian.
 *
 *  \return Constant GSL_SUCCESS to confirm that the computation was successful.
 */
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *ode_params)
{
  (void)(t);

  struct ODEParameters parameters = *(struct ODEParameters *)ode_params;
  double rho0                     = parameters.rho0;
  double g0                       = parameters.g0;
  double eta_ion                  = eval_interp(parameters.interp_ion, y[3]);
  double eta_diss                 = eval_interp(parameters.interp_diss, y[3]);

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 5, 5);
  gsl_matrix *m            = &dfdy_mat.matrix;

  gsl_matrix_set(m, 0, 0,
                 (0.03283531956573588 * rho0 * y[2] + 0.5 * rho0 * eta_ion * y[2] -
                  535200.887413877 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0]) /
                     (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(
      m, 0, 1, ((0.03283531956573588 + 0.5 * eta_ion) * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(
      m, 0, 2,
      (0.06567063913147177 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) + pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) * eta_ion +
       0.03283531956573588 * rho0 * y[2] + 0.5 * rho0 * eta_ion * y[2]) /
          (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(m, 0, 3, 0);
  gsl_matrix_set(m, 0, 4, 0);

  gsl_matrix_set(m, 1, 0,
                 (-0.5 * rho0 * eta_ion * y[2] + 0.5 * rho0 * eta_diss * y[2] +
                  535200.887413877 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0] -
                  0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] -
                  9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3]) /
                     (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(m, 1, 1,
                 (-0.5 * rho0 * eta_ion * y[2] + 0.5 * rho0 * eta_diss * y[2] -
                  0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0] -
                  0.2479695231261206 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] -
                  0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[2] -
                  9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0] * y[3] -
                  18505.188292994077 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3] -
                  9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[2] * y[3]) /
                     (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(m, 1, 2,
                 (pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) * eta_diss - 1.0 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) * eta_ion -
                  0.5 * rho0 * eta_ion * y[2] + 0.5 * rho0 * eta_diss * y[2] -
                  0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] -
                  9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3]) /
                     (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(m, 1, 3, (-11435.019367644041 * y[0] - 11435.019367644041 * y[1] - 11435.019367644041 * y[2]) * rho0 * y[1]);
  gsl_matrix_set(m, 1, 4, 0);

  gsl_matrix_set(
      m, 2, 0,
      (-0.5 * rho0 * y[2] - 0.5 * rho0 * eta_diss * y[2] + 0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] +
       9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3]) /
          (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(
      m, 2, 1,
      (-0.5 * rho0 * y[2] - 0.5 * rho0 * eta_diss * y[2] + 0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0] +
       0.2479695231261206 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] +
       0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[2] +
       9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[0] * y[3] +
       18505.188292994077 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3] +
       9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[2] * y[3]) /
          (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(
      m, 2, 2,
      (-1.0 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) - 1.0 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) * eta_diss -
       0.5 * rho0 * y[2] - 0.5 * rho0 * eta_diss * y[2] + 0.1239847615630603 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] +
       9252.594146497038 * sqrt((y[0] + y[1] + y[2]) * rho0) * rho0 * g0 * y[1] * y[3]) /
          (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(m, 2, 3, (11435.019367644041 * y[0] + 11435.019367644041 * y[1] + 11435.019367644041 * y[2]) * rho0 * y[1]);
  gsl_matrix_set(m, 2, 4, 0);

  gsl_matrix_set(m, 3, 0,
                 ((0.0037987669804047212 - 0.5 * y[3]) * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(m, 3, 1,
                 ((0.0037987669804047212 - 0.5 * y[3]) * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(
      m, 3, 2,
      (0.0075975339608094425 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) - 1.0 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) * y[3] +
       0.0037987669804047212 * rho0 * y[2] - 0.5 * rho0 * y[2] * y[3]) /
          (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(m, 3, 3, (-1 * sqrt((y[0] + y[1] + y[2]) * rho0) * y[2]) / (0.8091454722567165 * g0));
  gsl_matrix_set(m, 3, 4, 0);

  gsl_matrix_set(m, 4, 0, (0.4671646804342641 * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(m, 4, 1, (0.4671646804342641 * rho0 * y[2]) / (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(m, 4, 2,
                 (0.9343293608685282 * pow(sqrt((y[0] + y[1] + y[2]) * rho0), 2) + 0.4671646804342641 * rho0 * y[2]) /
                     (0.8091454722567165 * sqrt((y[0] + y[1] + y[2]) * rho0) * g0));
  gsl_matrix_set(m, 4, 3, 0);
  gsl_matrix_set(m, 4, 4, 0);

  dfdt[0] = 0;
  dfdt[1] = 0;
  dfdt[2] = 0;
  dfdt[3] = 0;
  dfdt[4] = 0;

  return GSL_SUCCESS;
};



/*! \brief Compute the stellar fraction according to our model.
 *
 *  This fraction can be used to calculate the probability of forming a star population.
 *
 *  \param[in] index Index of the gas cell in question.
 *  \param[in] dt Integration time in internal units.
 *
 *  \return The final stellar fraction ∈ [0, 1].
 */
double stellar_fraction(const int index, const double dt)
{
  double int_time = dt * T_GYR;

  /* Compute eta_ion and eta_diss with a interpolation function */
  glob_t globbuf;
  int glob_err = glob(IMF_PATH, GLOB_ERR, NULL, &globbuf);

  if(glob_err != 0)
    {
      printf("GLOB error: ");
      check_error(GLOB_ERROR, glob_err);
    }

  if(globbuf.gl_pathc < 2)
    {
      fprintf(stdout, "Glob could not find all eta files\n");
      fflush(stdout);
      endrun(8888);
    }

  struct InterpFunc interp_ion;
  struct InterpFunc interp_diss;

  interp_ion.acc    = gsl_interp_accel_alloc();
  interp_ion.interp = gsl_spline_alloc(gsl_interp_linear, globbuf.gl_pathc);
  interp_ion.x_min  = NAN;
  interp_ion.x_max  = NAN;
  get_eta(globbuf, int_time, 0, &interp_ion);

  interp_diss.acc    = gsl_interp_accel_alloc();
  interp_diss.interp = gsl_spline_alloc(gsl_interp_linear, globbuf.gl_pathc);
  interp_diss.x_min  = NAN;
  interp_diss.x_max  = NAN;
  get_eta(globbuf, int_time, 1, &interp_diss);

  struct ODEParameters ode_params;

  ode_params.g0          = 1.0;
  ode_params.interp_ion  = &interp_ion;
  ode_params.interp_diss = &interp_diss;
  ode_params.rho0 = 0.0;

  /* Total baryonic density in [Mₒ pc^(-3)] */
  double density = SphP[index].d.Density * RHO_COSMO;
  double ne = SphP[index].Ne, nh0 = 0, nHeII = 0;
  double ionized_fraction = 0.0, atomic_fraction = 0.0, molecular_fraction = 0.0, metal_fraction = 0.0, star_fraction = 0.0;
  double y[5], u, fac_entr_to_u;
  double T0 = 0.0;

  fac_entr_to_u = pow(SphP[index].d.Density * All.cf_a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
  //u             = DMAX(All.MinEgySpec, SphP[index].Entropy * fac_entr_to_u);
  //AbundanceRatios(u, SphP[index].d.Density * All.cf_a3inv, &ne, &nh0, &nHeII);
  
  /* Solve the ODE system using the BDF method */
  gsl_odeiv2_system sys     = {sf_ode, jacobian, 5, (void *)&ode_params};
  gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf, int_time * 1e-3, ABS_TOL, REL_TOL);
  int gsl_err;
  
  for(int i = 0; i < DIVISIONS; ++i)
    {
      y[0]        = 0.3; /* Ionized fraction */
      y[1]        = 0.3;     /* Atomic fraction */
      y[2]        = 0.4;     /* Molecular fraction */
      y[3]        = 1e-4;     /* Metal fraction */
      y[4]        = 0.0;     /* Stellar fraction */
	  
      ode_params.rho0 = density * F_RHO[i];
      gsl_odeiv2_driver_reset(driver);

      gsl_err = gsl_odeiv2_driver_apply(driver, &T0, int_time, y);
      if(gsl_err != GSL_SUCCESS)
        {
          printf("GSL ODE solver error: ");
          check_error(GLS_ODE_ERROR, gsl_err);
        }
	  
	  ionized_fraction += y[0] * PDF[i];
      atomic_fraction += y[1] * PDF[i];
      molecular_fraction += y[2] * PDF[i];
      metal_fraction += y[3] * PDF[i];
      star_fraction += y[4] * PDF[i];
	  T0 = 0.0;
    }
    
  gsl_odeiv2_driver_free(driver);
  gsl_spline_free(interp_ion.interp);
  gsl_spline_free(interp_diss.interp);
  gsl_interp_accel_free(interp_ion.acc);
  gsl_interp_accel_free(interp_diss.acc);

  /* Update internal variable storing the fractions */
  SphP[index].fHII  = ionized_fraction;
  SphP[index].fHI   = atomic_fraction;
  SphP[index].fmol  = molecular_fraction;
  SphP[index].fZ    = metal_fraction;
  SphP[index].fstar = star_fraction;

  /* Return the stellar mass fraction */
  return star_fraction;
}

//#endif /* #ifdef EZ_SFR */

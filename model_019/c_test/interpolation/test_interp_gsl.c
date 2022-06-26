/*
 * Test with:
 * gcc -O3 test_interp_gsl.c -IC:/msys64/mingw64/include -LC:/msys64/mingw64/lib -lgsl -lgslcblas -lm -o test_interp_gsl && ./test_interp_gsl.exe
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

double fun_0(double *idx)
{
    return idx[0];
}

double fun_1(double *idx)
{
    return idx[1];
}

double inv_fun_0(double *idx)
{
    return idx[0];
}

double inv_fun_1(double *idx)
{
    return idx[1];
}

void compareFiles(FILE *fp1, FILE *fp2)
{
    char ch1 = getc(fp1);
    char ch2 = getc(fp2);

    int error = 0, pos = 0, line = 1;

    while (ch1 != EOF && ch2 != EOF)
    {
        pos++;

        if (ch1 == '\n' && ch2 == '\n')
        {
            line++;
            pos = 0;
        }

        if (ch1 != ch2)
        {
            error++;
            printf("Line number : %d \tError position : %d \n", line, pos);
        }

        ch1 = getc(fp1);
        ch2 = getc(fp2);
    }

    printf("Total errors : %d\t", error);
}

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
static double **read_ftable(const char *filename, size_t *rows, size_t *cols)
{
    if (rows == NULL || cols == NULL || filename == NULL)
    {
        return NULL;
    }

    *rows = 0;
    *cols = 0;

    FILE *fp = fopen(filename, "r");

    if (fp == NULL)
    {
        printf("Could not open %s\n", filename);
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
            printf("The allocation of memory in read_ftable() failed\n");
            return matrix;
        }

        matrix = tmp;

        matrix[*rows] = (double *)calloc(*cols, sizeof(*matrix[*rows]));

        if (matrix[*rows] == NULL)
        {
            fclose(fp);
            printf("The allocation of memory in read_ftable() failed\n");
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

/*! \brief Evaluates a n-dimensional function using interpolation.
 *
 *  Return the value of a function F at a given point x, using F(x_i) for several known x_i.
 *
 *  \param[in] table Matrix with the known values
 *  \param[in] nrows Number of rows of `table`
 *  \param[in] ncols Number of columns of `table`
 *  \param[in] x Point at which the function is not known
 *  \param[in] ngrid Number of points in the grid for each dimension
 *  \param[in] fun List of function that transform indices to values in the grid
 *  \param[in] inv_fun List of function that transform values in the grid to indices
 *
 *  \return F evaluated at x.
 */
static double interpolate(double **table, size_t *nrows, size_t *ncols, double *x, int *ngrid, double (*fun[])(double *), double (*inv_fun[])(double *))
{
    const size_t DIMENSIONS = *ncols - 1;     // Number of dimensions of the problem
    const size_t N_VERT = pow(2, DIMENSIONS); // Number of vertices of a n-rectangle

    /***********************************************************************************************
     * Find and load the vertices of the n-rectangle around the target point
     **********************************************************************************************/

    double *idx_low = (double *)malloc(DIMENSIONS * sizeof(double));  // Index of the low vertex for each dimension
    double *idx_high = (double *)malloc(DIMENSIONS * sizeof(double)); // Index of the high vertex for each dimension
    double *val_low = (double *)malloc(DIMENSIONS * sizeof(double));  // Value of the low vertex for each dimension
    double *val_high = (double *)malloc(DIMENSIONS * sizeof(double)); // Value of the high vertex for each dimension

    double *max_vals = table[*nrows - 1]; // Maximum possible values for each dimension
    double *min_vals = table[0];          // Minimum possible values for each dimension

    for (size_t i = 0; i < DIMENSIONS; ++i)
    {
        if (x[i] < min_vals[i])
        {
            idx_low[i] = 0;
            idx_high[i] = 0;
        }
        else if (x[i] > max_vals[i])
        {
            idx_low[i] = ngrid[i] - 1;
            idx_high[i] = ngrid[i] - 1;
        }
        else
        {
            idx_low[i] = floor((*inv_fun[i])(x));
            idx_high[i] = ceil((*inv_fun[i])(x));
        }

        val_low[i] = (*fun[i])(idx_low);
        val_high[i] = (*fun[i])(idx_high);
    }

    double **vertices = (double **)malloc(N_VERT * sizeof(double *)); // Table entries for each vertex in the n-rectangle
    int *idx_vert = (int *)malloc(DIMENSIONS * sizeof(int));          // Local index for each dimension of a vertex
    int idx_global;                                                   // Global index of a vertex
    for (size_t i = 0; i < N_VERT; ++i)
    {
        vertices[i] = (double *)malloc(*ncols * sizeof(double));
        int l = i;
        int dim = 0;
        /* This runs up to the most significant 1 */
        while (l > 0)
        {
            if (l & 1)
            {
                /* If current bit is 1 */
                idx_vert[dim] = (int)idx_high[dim];
            }
            else
            {
                /* If current bit is 0 */
                idx_vert[dim] = (int)idx_low[dim];
            }

            l = l >> 1;
            dim++;
        }
        /* The rest of the bits are 0, so use idx_low */
        for (; dim < DIMENSIONS; ++dim)
        {
            idx_vert[dim] = (int)idx_low[dim];
        }

        /* Compute the global index of a vertex using the local ones */
        idx_global = idx_vert[DIMENSIONS - 1];
        for (size_t j = 0; j < DIMENSIONS - 1; ++j)
        {
            int mul_ngrid = 1;
            for (size_t k = j + 1; k < DIMENSIONS; ++k)
            {
                mul_ngrid *= ngrid[k];
            }
            idx_global += idx_vert[j] * mul_ngrid;
        }

        for (size_t j = 0; j <= DIMENSIONS; ++j)
        {
            vertices[i][j] = table[idx_global][j];
        }
    }

    free(idx_vert);
    free(idx_low);
    free(idx_high);

    /***********************************************************************************************
     * Compute the value for the target point using n-lineal interpolation
     * See Zhang et al. (2021), https://doi.org/10.1145/3423184 and references therein
     **********************************************************************************************/

    double interp_val = 0.0; // Value for the target point
    double weight;           // Interpolation weight (fractional volume of the sub n-rectangles)
    int degen_factor;        // Degeneracy factor (for when the target point lies in an edge or face of the n-rectangle)

    for (size_t i = 0; i < N_VERT; ++i)
    {
        weight = 1.0;
        degen_factor = 0;
        for (size_t j = 0; j < DIMENSIONS; ++j)
        {
            if (val_low[j] != val_high[j])
            {
                weight *= vertices[i][j] == val_low[j] ? (val_high[j] - x[j]) : (x[j] - val_low[j]);
                weight /= (val_high[j] - val_low[j]);
            }
            else
            {
                degen_factor++;
            }
        }

        interp_val += vertices[i][DIMENSIONS] * weight / pow(2, degen_factor);
    }

    for (size_t i = 0; i < N_VERT; ++i)
    {
        free(vertices[i]);
    }
    free(vertices);
    free(val_low);
    free(val_high);

    return interp_val;
}

int main()
{

    /***********************************************************************************************
     * Reference values
     ***********************************************************************************************/

    /* Set up GSL interpolation */
    const gsl_interp2d_type *T = gsl_interp2d_bilinear;
    const size_t N = 100;
    const double xa[] = {0.0, 1.0};
    const double ya[] = {0.0, 1.0};
    const size_t nx = sizeof(xa) / sizeof(double);
    const size_t ny = sizeof(ya) / sizeof(double);
    double *za = malloc(nx * ny * sizeof(double));
    gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
    gsl_interp_accel *xacc = gsl_interp_accel_alloc();
    gsl_interp_accel *yacc = gsl_interp_accel_alloc();
    size_t i, j;

    /* Set z grid values */
    gsl_spline2d_set(spline, za, 0, 0, 0.0);
    gsl_spline2d_set(spline, za, 0, 1, 1.0);
    gsl_spline2d_set(spline, za, 1, 0, 1.0);
    gsl_spline2d_set(spline, za, 1, 1, 0.5);
    /* Write values to a file */
    FILE *fp;
    fp = fopen("./data.txt", "w+");
    fprintf(fp, "0 0 0.0\n0 1 1.0\n1 0 1.0\n1 1 0.5\n");
    fclose(fp);

    /* Initialize interpolation */
    gsl_spline2d_init(spline, xa, ya, za, nx, ny);

    /* Interpolate N values in x and y and write the result to a file */
    FILE *fp_gsl;
    fp_gsl = fopen("./gsl_out.txt", "w+");
    for (i = 0; i < N; ++i)
    {
        double xi = i / (N - 1.0);

        for (j = 0; j < N; ++j)
        {
            double yj = j / (N - 1.0);
            double zij = gsl_spline2d_eval(spline, xi, yj, xacc, yacc);

            fprintf(fp_gsl, "%f %f %f\n", xi, yj, zij);
        }
        fprintf(fp_gsl, "\n");
    }

    gsl_spline2d_free(spline);
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);
    free(za);

    /***********************************************************************************************
     * Multivariate interpolation
     ***********************************************************************************************/

    /* Load known data */
    size_t cols, rows;
    double **data_table;
    char paths[] = "./data.txt";
    data_table = read_ftable(paths, &rows, &cols);

    /* Set up interpolation */
    int n[] = {2, 2};
    double (*fun[])(double *) = {fun_0, fun_1};
    double (*inv_fun[])(double *) = {inv_fun_0, inv_fun_1};

    /* Interpolate N values in x and y and write the result to a file */
    FILE *fp_result;
    fp_result = fopen("./mi_out.txt", "w+");
    for (i = 0; i < N; ++i)
    {
        double xi = i / (N - 1.0);

        for (j = 0; j < N; ++j)
        {
            double yj = j / (N - 1.0);
            double point[] = {xi, yj};
            double zij = interpolate(data_table, &rows, &cols, point, n, fun, inv_fun);
            fprintf(fp_result, "%f %f %f\n", xi, yj, zij);
        }
        fprintf(fp_result, "\n");
    }

    compareFiles(fp_gsl, fp_result);

    fclose(fp_gsl);
    fclose(fp_result);

    return 0;
}
/*
 * Make into a shared library in Windows with: gcc -fpic -shared -O3 interpolation.c -o libinterpolation.dll
 * Make into a shared library in Linux with: gcc -fpic -shared -O3 interpolation.c -o libinterpolation.so
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PRECISION 1.0e-10
#define DEQU(a, b) (fabs(a - b) < PRECISION)
#define DNEQ(a, b) (fabs(a - b) >= PRECISION)

/* Index to value functions */
static double fun_if(double *idx) { return (0.01 + 0.05157894736842105 * idx[0]); }
static double fun_zf(double *idx) { return (0.002105263157894737 * idx[1]); }
static double fun_rh(double *idx) { return (0.002471402662718477 + 0.001821033540950457 * idx[2]); }
static double fun_it(double *idx) { return (0.001 + 0.0004736842105263158 * idx[3]); }
/* Value to index functions */
static double inv_fun_if(double *idx) { return (-0.19387755102040816 + 19.387755102040817 * idx[0]); }
static double inv_fun_zf(double *idx) { return (475.0 * idx[1]); }
static double inv_fun_rh(double *idx) { return (-1.357142857142857 + 549.138704758874 * idx[2]); }
static double inv_fun_it(double *idx) { return (-2.1111111111111107 + 2111.111111111111 * idx[3]); }

double (*FUN[])(double *) = {fun_if, fun_zf, fun_rh, fun_it};
double (*INV_FUN[])(double *) = {inv_fun_if, inv_fun_zf, inv_fun_rh, inv_fun_it};
const int NGRID[] = {20, 20, 20, 20};

/*! \brief This function loads a text file into memory as a C array of arrays.
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
 *
 *  \return F evaluated at x.
 */
double interpolate(double **table, size_t *nrows, size_t *ncols, double *x)
{
    const size_t DIMENSIONS = *ncols - 1;     // Number of dimensions of the problem
    const size_t N_VERT = pow(2, DIMENSIONS); // Number of vertices of a n-rectangle

    /*************************************************************************************************
     * Find and load the vertices of the n-rectangle around the target point
     ************************************************************************************************/

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
            idx_low[i] = NGRID[i] - 1;
            idx_high[i] = NGRID[i] - 1;
        }
        else
        {
            idx_low[i] = floor((*INV_FUN[i])(x));
            idx_high[i] = ceil((*INV_FUN[i])(x));
        }

        val_low[i] = (*FUN[i])(idx_low);
        val_high[i] = (*FUN[i])(idx_high);

        /* If the point coincides with a vertex, collapse that dimension */
        if (DEQU(x[i], val_low[i]) && DNEQ(x[i], val_high[i]))
        {
            val_high[i] = val_low[i];
            idx_high[i] = idx_low[i];
        }
        else if (DEQU(x[i], val_high[i]) && DNEQ(x[i], val_low[i]))
        {
            val_low[i] = val_high[i];
            idx_low[i] = idx_high[i];
        }
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
                mul_ngrid *= NGRID[k];
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

    /*************************************************************************************************
     * Compute the value for the target point using n-lineal interpolation
     * See Zhang et al. (2021), https://doi.org/10.1145/3423184 and references therein
     ************************************************************************************************/

    double interp_val = 0.0; // Value for the target point
    double weight;           // Interpolation weight (fractional volume of the sub n-rectangles)
    int degen_factor;        // Degeneracy factor (for when the target point lies in an edge or face of the n-rectangle)

    for (size_t i = 0; i < N_VERT; ++i)
    {
        weight = 1.0;
        degen_factor = 0;
        for (size_t j = 0; j < DIMENSIONS; ++j)
        {
            if (DNEQ(val_low[j], val_high[j]))
            {
                weight *= DEQU(vertices[i][j], val_low[j]) ? (val_high[j] - x[j]) : (x[j] - val_low[j]);
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
/*
 * Make into a shared library in Windows with: gcc -fpic -shared -O3 -march=native -mtune=native interpolation_private.c -o libinterpolationprivate.dll
 * Make into a shared library in Linux with: gcc -fpic -shared -O3 -march=native -mtune=native interpolation_private.c -o libinterpolationprivate.so
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEQU(a, b) (fabs(a - b) < 10.0e-11)
#define DNEQ(a, b) (fabs(a - b) >= 10.0e-11)

#define N_COLS 5       // Number of columns in the table
#define N_ROWS 1350000 // Number of rows in the table
#define N_DIMS 4       // Number of dimensions of the problem
#define N_VERT 16      // Number of vertices of a n-rectangle = 2^N_DIMS

/* Index to value functions */
static double fun_if(double *idx) { return (0.2 + 0.027586206896551727 * idx[0]); }
static double fun_zf(double *idx) { return (exp(-11.22025585100741 + 0.3314075151232176 * idx[1])); }
static double fun_rh(double *idx) { return (exp(-4.961845129926824 + 0.260936553893533 * idx[2])); }
static double fun_it(double *idx) { return (1.0e-5 + 3.103448275862069e-6 * idx[3]); }

/* Value to index functions */
static double inv_fun_if(double *idx) { return (-7.250000000000001 + 36.25 * idx[0]); }
static double inv_fun_zf(double *idx) { return (33.8563712015876 + 6.947896435414094 * log10(idx[1])); }
static double inv_fun_rh(double *idx) { return (19.015523336570737 + 8.824310195855288 * log10(idx[2])); }
static double inv_fun_it(double *idx) { return (-3.2222222222222223 + 322222.22222222225 * idx[3]); }

static double (*FUN[])(double *) = {fun_if, fun_zf, fun_rh, fun_it};
static double (*INV_FUN[])(double *) = {inv_fun_if, inv_fun_zf, inv_fun_rh, inv_fun_it};

static const int NGRID[] = {30, 30, 50, 30};

/*! \brief Loads a text file into memory as a C array.
 *
 *  Read a text file with space separated values, and stores the numerical
 *  data as doubles in a C array. Each line in the file must be at most 1024
 *  characters long.
 *
 *  \param[in] filepath Path to the file.
 *
 *  \return A pointer to the array storing the data, or NULL if the request failed.
 */
double *read_ftable(const char *filepath)
{
    FILE *file_ptr = fopen(filepath, "r");
    if (file_ptr == NULL)
    {
        printf("Could not open %s\n", filepath);
        return NULL;
    }

    double *table = malloc(N_ROWS * N_COLS * sizeof(double));
    if (table == NULL)
    {
        printf("The allocation of memory in read_ftable() failed\n");
        fclose(file_ptr);
        return NULL;
    }

    char line[1024];
    for (size_t i = 0; i < N_ROWS; ++i)
    {
        if (fgets(line, 1024, file_ptr) != NULL)
        {
            char *scan = line;
            int offset = 0;
            for (size_t j = 0; j < N_COLS; ++j)
            {
                if (sscanf(scan, "%lf%n", &(table[i * N_COLS + j]), &offset) == 1)
                {
                    scan += offset;
                }
                else
                {
                    printf("Could not read some character in %s\n", filepath);
                    fclose(file_ptr);
                    return NULL;
                }
            }
        }
        else
        {
            printf("Could not read some line in %s\n", filepath);
            fclose(file_ptr);
            return NULL;
        }
    }

    fclose(file_ptr);
    return table;
}

/*! \brief Evaluates a n-dimensional function using interpolation.
 *
 *  Return the value of a function F at a given point x, using known F(x_i) for several x_i.
 *
 *  \param[in] table Matrix with the known values
 *  \param[in] x Point at which the function is not known
 *
 *  \return F evaluated at x.
 */
double interpolate(double *table, double *x)
{
    /*************************************************************************************************
     * Find the vertices of the n-rectangle around the target point
     ************************************************************************************************/

    double idx_low[N_DIMS];  // Index of the low vertex for each dimension
    double idx_high[N_DIMS]; // Index of the high vertex for each dimension
    double val_low[N_DIMS];  // Value of the low vertex for each dimension
    double val_high[N_DIMS]; // Value of the high vertex for each dimension

    for (size_t i = 0; i < N_DIMS; ++i)
    {
        if (x[i] < table[i])
        {
            idx_low[i] = 0;
            idx_high[i] = 0;
        }
        else if (x[i] > table[(N_ROWS - 1) * N_COLS + i])
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

    /*************************************************************************************************
     * Find the rows corresponding to the vertices
     ************************************************************************************************/

    double vertices[N_VERT * N_COLS]; // Entries in the table for each vertex of the n-rectangle
    size_t idx_vert[N_DIMS];          // Local index for each dimension of a vertex
    size_t idx_global;                // Global index of a vertex
    for (size_t i = 0; i < N_VERT; ++i)
    {
        size_t bit_i = i;
        size_t dimension = 0;

        /* This runs up to the most significant 1 */
        while (bit_i > 0)
        {
            if (bit_i & 1)
            {
                /* If current bit is 1 */
                idx_vert[dimension] = (size_t)idx_high[dimension];
            }
            else
            {
                /* If current bit is 0 */
                idx_vert[dimension] = (size_t)idx_low[dimension];
            }

            bit_i = bit_i >> 1;
            dimension++;
        }

        /* The rest of the bits are 0, so use idx_low */
        for (; dimension < N_DIMS; ++dimension)
        {
            idx_vert[dimension] = (size_t)idx_low[dimension];
        }

        /* Compute the global index of a vertex */
        idx_global = idx_vert[N_DIMS - 1];
        for (size_t j = 0; j < N_DIMS - 1; ++j)
        {
            size_t mul_ngrid = 1;
            for (size_t k = j + 1; k < N_DIMS; ++k)
            {
                mul_ngrid *= NGRID[k];
            }
            idx_global += idx_vert[j] * mul_ngrid;
        }

        for (size_t j = 0; j < N_COLS; ++j)
        {
            vertices[i * N_COLS + j] = table[idx_global * N_COLS + j];
        }
    }

    /*************************************************************************************************
     * Compute the value for the target point using n-linear interpolation
     * See Zhang et al. (2021), https://doi.org/10.1145/3423184 and references therein
     ************************************************************************************************/

    double interp_val = 0.0;   // Value for the target point
    double weight = 1.0;       // Interpolation weight (fractional volume of the sub n-rectangles)
    size_t degen_factor = 0.0; // Degeneracy factor (for when the target point lies in an edge or face of the n-rectangle)

    for (size_t i = 0; i < N_VERT; ++i)
    {
        for (size_t j = 0; j < N_DIMS; ++j)
        {
            if (DNEQ(val_low[j], val_high[j]))
            {
                weight *= DEQU(vertices[i * N_COLS + j], val_low[j]) ? (val_high[j] - x[j]) : (x[j] - val_low[j]);
                weight /= (val_high[j] - val_low[j]);
            }
            else
            {
                degen_factor++;
            }
        }

        interp_val += vertices[i * N_COLS + N_DIMS] * weight / (0 | 1UL << degen_factor);
        weight = 1.0;
        degen_factor = 0;
    }

    return interp_val;
}

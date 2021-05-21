#include "libtsnnls-2.3.4/tsnnls/lsqr.h"
#include "libtsnnls-2.3.4/tsnnls/tsnnls.h"
#include "tsnnls_wrappers.h"
#include <stdio.h>

/*Creates taucs matrix representation for sparse matrix*/
taucs_ccs_matrix* create_matrix(int* rows, int* cols, double* vals, double* b, int num_rows, int num_cols, int nnz)
{
	double *A_vals = (double *)(calloc(num_rows*num_cols,sizeof(double)));

	for(unsigned int i = 0; i < nnz; i++)
	{
		A_vals[rows[i]*num_cols+cols[i]] = vals[i];
	}

	taucs_ccs_matrix *A = taucs_construct_sorted_ccs_matrix(A_vals, num_cols, num_rows);
	
	free(A_vals);

	return A;
}

/*Uses tsnnls to solve a non-negative least-squares problem*/
double* nnls_wrapper(int* rows, int* cols, double* vals, double* b, int num_rows, int num_cols, int nnz, double tol)
{
	double residual;
	
	taucs_ccs_matrix* A = create_matrix(rows, cols, vals, b, num_rows, num_cols, nnz);
	
	double* x = t_snnls(A, b, &residual, tol, 1);
	
	free(A);

	return x;
}

/*Uses tsnnls to solve a standard least-squares problem*/
double* lsqr_wrapper(int* rows, int* cols, double* vals, double* b, int num_rows, int num_cols, int nnz, double tol)
{
	taucs_ccs_matrix* A = create_matrix(rows, cols, vals, b, num_rows, num_cols, nnz);
	double* x = t_lsqr(A, b);
	free(A);

	return x;
}

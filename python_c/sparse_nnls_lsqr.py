import scipy
import numpy as np
from scipy.sparse import *
from ctypes import *
from numpy.ctypeslib import ndpointer
from scipy.optimize import *

#tsnnls: http://www.jasoncantarella.com/wordpress/software/tsnnls/
#compilation: cc -fPIC -shared -o tsnnls_wrappers.so tsnnls_wrappers.c tsnnls_wrappers.h -ltsnnls

so_file = "/home/arlei/DTRA/code/tsnnls_wrappers.so"

def extract_matrix_nnz(A):
	'''
		Getting non-zeros from the sparse matrix A.
	'''
	num_rows = A.shape[0]
	num_cols = A.shape[1]
	nnz = A.nnz
	
	rows, cols = A.nonzero()
	rows = list(rows)
	cols = list(cols)
	vals = []

	for i in range(nnz):
		r = rows[i]
		c = cols[i]
		vals.append(A[r,c])
	
	return num_rows, num_cols, nnz, rows, cols, vals

def sparse_nnls(A, b, tol=1e-10, eps=1e-5):
	'''
		Solving non-negative least squares by
		using the tsnnls library (in C).
	'''
	num_rows, num_cols, nnz, rows, cols, vals = extract_matrix_nnz(A)
	A_sum = np.array(A.sum(axis=1).T)[0]
	c_b = (c_double * b.shape[0])(*list(b+eps*A_sum))
	c_rows = (c_int * nnz)(*rows)
	c_cols = (c_int * nnz)(*cols)
	c_vals = (c_double * nnz)(*vals)
	c_tol = (c_double) (tol)

	tsnnls_wrappers = CDLL(so_file)
	tsnnls_wrappers.nnls_wrapper.restype = ndpointer(dtype=c_double, shape=(num_cols))
	
	max_iter = 20
	i = 0
	c_x = None
	while c_x is None and i < max_iter:
		#This doesnt look nice, but seem to work
		#We have modified the tsnnl implementation to return
		#NULL instead of an assert in a particular scenario
		#and adding a non-negative vector to b seems to avoid 
		#that while adding some error to the result.
		#We try to add as little error as possible.
		try:
			c_x = tsnnls_wrappers.nnls_wrapper(c_rows, c_cols, c_vals, c_b, num_rows, num_cols, nnz, c_tol)
			neg = c_x <= eps
			c_x[neg] = 0
		except ValueError:
			#NULL
			eps = eps * 2
			c_b = (c_double * b.shape[0])(*list(b+eps*A_sum))
			i = i + 1
	
	if c_x is None:
		print("Warning: TNNLS error.")
		return np.zeros(num_cols)
	else:
		return c_x

def sparse_lsqr(A, b):
	'''
		Solving standard least squares by
		using the tsnnls library (in C).
	'''
	num_rows, num_cols, nnz, rows, cols, vals = extract_matrix_nnz(A)

	c_b = (c_double * b.shape[0])(*list(b))
	c_rows = (c_int * nnz)(*rows)
	c_cols = (c_int * nnz)(*cols)
	c_vals = (c_double * nnz)(*vals)

	tsnnls_wrappers = CDLL(so_file)
	tsnnls_wrappers.lsqr_wrapper.restype = ndpointer(dtype=c_double, shape=(num_cols))
	
	try:
		#tsnnls results seem to flip the sign for some reason
		c_x = tsnnls_wrappers.lsqr_wrapper(c_rows, c_cols, c_vals, c_b, num_rows, num_cols, nnz)
	except ValueError:
		print("TNNLS error.")
		c_x = None
	
	return c_x

#A = scipy.sparse.load_npz("A.npz")
#b = np.load("b.npz.npy")

#print(A.shape)

#new_A = A[A.getnnz(1) > 0]
#new_b = b[A.getnnz(1) > 0]

#print(new_A.shape)

#for j in range(new_A.shape[1]):
#	if new_A[:,j].sum() == 0:
#		print(j)

#print(np.linalg.matrix_rank(A.todense()))
#print(A.shape)

#A = np.random.random((5,5))
#b = np.random.random(5)

#print(np.linalg.cond(A.todense()))

#print(sparse_nnls(scipy.sparse.csr_matrix(new_A), -new_b, 100))

#A = np.random.random((5,5))
#b = np.random.random(5)
#print(sparse_lsqr(scipy.sparse.csr_matrix(A), -b))

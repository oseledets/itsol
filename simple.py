from _itsol_cffi import lib, ffi
import numpy as np
import scipy.sparse as sp

def get_sp_from_itsol(L_csr, mat):
    N = L_csr.n
    nnz = L_csr.indptr[N]
    int_size = mat.indices.dtype.itemsize
    double_size = mat.data.dtype.itemsize
    L_data = np.frombuffer(ffi.buffer(L_csr.data, nnz*double_size), dtype=np.float64)
    L_indices = np.frombuffer(ffi.buffer(L_csr.indices, nnz*int_size), dtype=np.int32)
    L_indptr = np.frombuffer(ffi.buffer(L_csr.indptr, (N+1)*int_size), dtype=np.int32)
    return L_data, L_indices, L_indptr

#Convert L to csr format, need to get back from buffer
    


sparse_mat = ffi.new("ITS_SparMat *")
sparse_lu = ffi.new("ITS_ILUSpar *")
n = 128
e = sp.eye(n)
z = sp.diags(np.ones(n-1), -1)
lp = 2*e - z - z.T
nnz_count = lp.indptr[1:]-lp.indptr[:-1]
lp2 = sp.kron(lp, e) + sp.kron(e, lp)
mat = lp2

#sparse_mat.nzcount = ffi.from_buffer("int[]", nnz_count)
#
#sparse_mat.n = n
#sparse_mat.ja = ffi.new(f"int*[{n}]")
#sparse_mat.ma = ffi.new(f"double*[{n}]")
#for i in range(n):
#    cur_indices = lp.indices[lp.indptr[i]:lp.indptr[i+1]]
#    cur_data = lp.data[lp.indptr[i]:lp.indptr[i+1]]
#    sparse_mat.ja[i] = ffi.from_buffer("int[]", cur_indices)
 #   sparse_mat.ma[i] = ffi.from_buffer("double[]", cur_data)



sparse_lu = lib.itsol_from_numpy(mat.shape[0], 2, ffi.from_buffer("double[]", mat.data), ffi.from_buffer("int[]", mat.indices), ffi.from_buffer("int[]", mat.indptr))

L_csr = lib.itsol_SpaFmtNumpy(sparse_lu.L) 
U_csr = lib.itsol_SpaFmtNumpy(sparse_lu.U)

L_data, L_indices, L_indptr = get_sp_from_itsol(L_csr, mat)
U_data, U_indices, U_indptr = get_sp_from_itsol(U_csr, mat)
U = sp.csr_matrix((U_data, U_indices, U_indptr), shape=[mat.shape[0], mat.shape[1]])
L = sp.csr_matrix((L_data, L_indices, L_indptr), shape=[mat.shape[0], mat.shape[1]])


#Now test prec solve as given by itsol
#y -- input, x -- output
rhs = np.ones(mat.shape[0])
sol = np.zeros(mat.shape[0])
lib.itsol_lusolC(ffi.from_buffer("double[]", rhs), ffi.from_buffer("double[]", sol), sparse_lu)

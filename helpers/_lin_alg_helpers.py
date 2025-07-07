def gen_invert_matrix(self, use_thread_group=False):
    """
    This function generates a matrix inversion function for cuda.
    The function employs Gaussian elimination.
    """

    self.gen_add_func_doc("Compute the inverse of a matrix", ["Uses gaussian elimination"], \
                          ['dimA is number of rows in A', \
                           'A is a pointer to the original invertible matrix. It is turned into an identity matrix', \
                           'Ainv is a pointer to an identity matrix that will be transformed into the inverse of A', \
                            's_temp is a pointer to temporary memory of size 4*dimA'])
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void invert_matrix(uint32_t dimA, T *A, T *Ainv, T *s_temp) {", True)
    self.gen_add_code_line("for (unsigned pivRC = 0; pivRC < dimA; pivRC++) {", True)   # iterate over diagonal
    self.gen_add_code_line("unsigned pivColOffset = pivRC*dimA;")                       # offset to first value of pivot row
    self.gen_add_code_line("T pvInv = static_cast<T>(1)/A[pivRC + pivColOffset];")      # 1/pivot

    # save the pivot row and column values
    self.gen_add_parallel_loop("ind", "dimA", use_thread_group)
    self.gen_add_code_line("s_temp[ind] = static_cast<T>(A[pivRC * dimA + ind]);")
    self.gen_add_code_line("s_temp[ind+dimA] = static_cast<T>(Ainv[pivRC * dimA + ind]);")
    self.gen_add_code_line("s_temp[ind+dimA*2] = static_cast<T>(A[pivRC + dimA * ind]);")
    self.gen_add_code_line("s_temp[ind+dimA*3] = static_cast<T>(Ainv[pivRC + dimA * ind]);")
    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)

    # run gaussian elimation for the pivot row and column
    self.gen_add_parallel_loop("ind", "dimA*dimA", use_thread_group)
    self.gen_add_code_line("unsigned row = ind % dimA, col = ind / dimA;")
    # apply to the pivot row
    self.gen_add_code_line("if (row == pivRC) {", True)
    self.gen_add_code_line("A[row * dimA + col] *= pvInv;") # put 1 on the diagonal by multiplying row by inverse
    self.gen_add_code_line("Ainv[row * dimA + col] *= pvInv;")
    self.gen_add_end_control_flow()
    # apply to other rows by reducing entries on the pivot column to 0s
    self.gen_add_code_line("else {", True)
    self.gen_add_code_line("T multiplier = s_temp[row+dimA*2] / s_temp[pivRC];")
    self.gen_add_code_line("A[row * dimA + col] -= multiplier * s_temp[col];")
    self.gen_add_code_line("Ainv[row * dimA + col] -= multiplier * s_temp[col+dimA];")
    self.gen_add_end_control_flow()

    self.gen_add_end_control_flow()
    self.gen_add_sync(use_thread_group)
    self.gen_add_end_control_flow()
    self.gen_add_end_function()
    return


def gen_matmul(self):
    """
    Generates the matrix multiplication helper function.
    This function allows for a transpose of B
    """
    self.gen_add_func_doc("Matrix multiplication helper function of AB", [], \
                          ['index - the index of the result vector', \
                           'A - pointer to the first matrix', \
                           'B - pointer to the second matrix', \
                           'dest - pointer to the destination matrix', \
                           'num - 36 or 6 depending on the indexing scheme', \
                           't - true => multiply with the transpose of B'])
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void matmul(int index, T *A, T *B, T *dest, int num, bool t) {", True)
    self.gen_add_code_line("int cur = 36*((index/num)%NUM_JOINTS);")
    self.gen_add_code_line("T *vec1 = &B[cur + (t*5+1)*(index%6)];")
    self.gen_add_code_line("T *vec2 = &A[6*(index/6)];")
    self.gen_add_code_line("dest[index] = dot_prod<T,6, 6, 1>(vec1, vec2);")
    self.gen_add_end_function()


def gen_matmul_trans(self):
    """
    Generates the matrix multiplication helper function where one of the
    matrices is transposed. Both A and B are 6x6 matrices.
    """
    self.gen_add_func_doc("Matrix multiplication helper function where one of the matrices is tranposed.", [], \
                          ['index - the index of the result vector', \
                           'A - pointer to the first 6x6 matrix', \
                           'B - pointer to the second 6x6 matrix', \
                           'dest - pointer to the destination matrix', \
                           'char trans_mat - a for A^TB, b for AB^T'])
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void matmul_trans(int index, T *A, T *B, T *dest, char trans_mat) {", True)
    self.gen_add_code_line("T *vec1;")
    self.gen_add_code_line("T *vec2;")
    self.gen_add_code_line("if (trans_mat == 'a'){", True)
    self.gen_add_code_line("vec1 = &A[6*(index%6)];")
    self.gen_add_code_line("vec2 = &B[6*(index/6)];")
    self.gen_add_code_line("dest[index] = dot_prod<T,6,1,1>(vec1, vec2);")
    self.gen_add_end_control_flow()
    self.gen_add_code_line("if (trans_mat == 'b'){", True)
    self.gen_add_code_line("vec1 = &A[index%6];")
    self.gen_add_code_line("vec2 = &B[index/6];")
    self.gen_add_code_line("dest[index] = dot_prod<T,6,6,6>(vec1, vec2);")
    self.gen_add_end_control_flow()
    self.gen_add_end_function()


def gen_outer_product(self):
    """
    This function generates the cuda for the outerProduct
    function.
    """
    self.gen_add_func_doc("Compute the outer product between two vectors: dest = ab^T", \
                          ["Function assumes it is called by a single thread."], \
                          ['a - first vector', \
                           'b - second vector', \
                           'dest - destination matrix', \
                           'aLength - length of a', \
                           'bLength - length of b', \
                           'idx - index of resulting matrix to be computed by this thread'])
    self.gen_add_code_line("template <typename T>")
    self.gen_add_code_line("__device__")
    self.gen_add_code_line("void outerProduct(T *a, T *b, T *dest, int aLength, int bLength, int idx) {", True)
    self.gen_add_code_line("int row = idx / bLength;")
    self.gen_add_code_line("int col = idx % bLength;")
    self.gen_add_code_line("if (row < aLength && col < bLength) dest[col * aLength + row] = a[row] * b[col];")
    self.gen_add_end_function()

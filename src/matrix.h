/////////////////////////////////////////////////////////////////////////////////////////
#ifndef __MATRIX_H__
#define __MATRIX_H__
/////////////////////////////////////////////////////////////////////////////////////////

// Compute an integral generator of the kernel of A
std::vector<std::vector<int>>
fourier_motzkin_kernel(const std::vector<std::vector<int>>& A);

// read a matrix from file (dense format)
bool 
read_mat(const char* fname,
         std::vector<std::vector<int>>& A);

// read a matrix from file, in GreatSPN sparse format
bool
read_pin(const char* fname, const size_t m,
         std::vector<std::vector<int>>& A);

std::vector<std::vector<int>>
transpose(const std::vector<std::vector<int>>& A);

std::vector<std::vector<int>>
make_symmetric(const std::vector<std::vector<int>>& A);

// for each column j<m find the first row i s.t. A[0..i-1, j] <= 0, and A[i,j] > 0
std::vector<size_t>
step_for_negative_removal(const std::vector<std::vector<int>>& A);

// move A in row echelon form
void row_echelon_form(std::vector<std::vector<int>>& A, 
                      std::vector<size_t> *leading_cols=nullptr);

void row_footprint_form(std::vector<std::vector<int>>& A, 
                        std::vector<size_t> *leading_cols=nullptr,
                        std::vector<size_t> *trailing_cols=nullptr);

size_t irank(const std::vector<std::vector<int>>& A);

// generating set of the nullspace of A
std::vector<std::vector<int>>
nullspace(std::vector<std::vector<int>>& A);

// generate the nullspace starting from the hnf/echelon form
std::vector<std::vector<int>>
nullspace_from_hnf(const std::vector<std::vector<int>>& A,
                   const std::vector<size_t>& pivot_cols);

// Count the space occupied by each column in a compact representation
void account_col_space(std::vector<size_t>& col_spaces, const int *const data);

// print matrix
void print_mat(const std::vector<std::vector<int>>& A,
               bool highlight_spans = false);

bool
integral_kernel_old(const std::vector<std::vector<int>>& A,
                    std::vector<std::vector<int>>& result,
                    bool verbose = false);

void 
hermite_normal_form(const std::vector<std::vector<int>>& A, 
                    std::vector<std::vector<int>>& H,
                    std::vector<std::vector<int>>& U,
                    std::vector<size_t> *leading_cols,
                    bool verbose = false);

// compute a set of Z-generators for the lattice of the kernel of A
std::vector<std::vector<int>>
integral_kernel_Zgens(const std::vector<std::vector<int>>& A,  bool verbose);

/////////////////////////////////////////////////////////////////////////////////////////
#endif // __MATRIX_H__
/////////////////////////////////////////////////////////////////////////////////////////

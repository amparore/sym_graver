#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
#include <cassert>
#include <algorithm>
#include <numeric>

#include "math_utils.h"
#include "matrix.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////

overflow_exception::overflow_exception(const char* message)
: message(message) {}

overflow_exception::overflow_exception(const std::string& message)
: message(message) {}

const char* overflow_exception::what() const throw() 
{ return message.c_str(); }

/////////////////////////////////////////////////////////////////////////////////////////

inline size_t count_nonzero(const std::vector<int>& v) {
    return std::count_if(v.begin(), v.end(), [](int x){ return x != 0;});
}

// return the index of the first non-zero entry in v
inline size_t leading(const std::vector<int>& v) {
    auto it = std::find_if(v.begin(), v.end(), [](int x){ return x != 0;});
    return (it==v.end()) ? size_t(-1) : it-v.begin();
}

// return the index of the last non-zero entry in v
inline size_t trailing(const std::vector<int>& v) {
    auto it = std::find_if(v.rbegin(), v.rend(), [](int x){ return x != 0;});
    return (it==v.rend()) ? size_t(-1) : v.size() - 1 - (it - v.rbegin());
}


// ensure that gcd(v) is 1
int reduce_gcd(std::vector<int>& v) {
    int g = 0;
    for (size_t i=0; i<v.size(); i++)
        if (v[i] != 0)
            g = gcd(g, abs(v[i]));
    if (g != 0) {
        for (size_t i=0; i<v.size(); i++)
            if (v[i] != 0)
                v[i] /= g;
    }
    return g;
}

// ensure that gcd of column j is 1
int reduce_column_gcd(std::vector<std::vector<int>>& A, size_t j) {
    int g = 0;
    for (size_t i=0; i<A.size(); i++)
        if (A[i][j] != 0)
            g = gcd(g, abs(A[i][j]));
    if (g != 0) {
        for (size_t i=0; i<A.size(); i++)
            if (A[i][j] != 0)
                A[i][j] /= g;
    }
    return g;
}

// // ensure that gcd(v) is 1 and sign is canonical
// void canonicalize(std::vector<int>& v) {
//     size_t l = leading(v);
//     if (l != size_t(-1)) {
//         if (v[l] < 0) { // change sign
//             for (size_t i=l; i<v.size(); i++) 
//                 v[i] = -v[i];
//         }
//         reduce_gcd(v);
//     }
// }

/////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////

// Count the space occupied by each column in a compact representation
void account_col_space(std::vector<size_t>& col_spaces, const int *const data) {
    for (size_t i=0; i<col_spaces.size(); i++) {
        size_t sp = 1;
        if (data[i] < 0) sp++;
        if (abs(data[i]) >= 10) sp++;
        if (abs(data[i]) >= 100) sp++;
        if (abs(data[i]) >= 1000) sp++;
        col_spaces[i] = max(col_spaces[i], sp);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

void print_mat(const std::vector<std::vector<int>>& A,
               bool highlight_spans) 
{
    if (A.empty())
        return;
    const size_t m = A.front().size();
    std::vector<size_t> col_spaces(m);
    for (const auto& row : A) {
        account_col_space(col_spaces, row.data());
    }
    for (const auto& row : A) {
        size_t l = leading(row), t = trailing(row);
        for (size_t i=0; i<row.size(); i++) {
            if (highlight_spans && l!=size_t(-1) && (i<l || i>t))
                cout << setw(col_spaces[i]) << "." << " ";
            else
                cout << setw(col_spaces[i]) << row[i] << " ";
        }
        cout << endl;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

// struct row_t {
//     std::vector<int> vd, vc;

//     // empty initial row
//     inline row_t(size_t n, size_t m) : vd(n), vc(m) { }
//     // sum of two rows
//     inline row_t(int mult1, const row_t& row1, int mult2, const row_t& row2);

//     void add_mult(int mult, const row_t& row);
//     void add_vec(const std::vector<int>& vec);

//     int gcd_dc();
//     void canonicalize();
// };

// row_t::row_t(int mult1, const row_t& row1, int mult2, const row_t& row2) 
// /**/ : row_t::row_t(row1.vd.size(), row1.vc.size()) 
// {
//     add_mult(mult1, row1);
//     add_mult(mult2, row2);
// }

// void row_t::add_mult(int mult, const row_t& row) {
//     assert(vd.size()==row.vd.size() && vc.size()==row.vc.size());
//     for (int j=0; j<vd.size(); j++)
//         vd[j] = add_exact(vd[j], multiply_exact(mult, row.vd[j]));
//     for (int j=0; j<vc.size(); j++)
//         vc[j] = add_exact(vc[j], multiply_exact(mult, row.vc[j]));
// }

// void row_t::add_vec(const std::vector<int>& vec) {
//     assert(vc.size()==vec.size());
//     for (int j=0; j<vc.size(); j++)
//         vc[j] = add_exact(vc[j], vec[j]);
// }

// int row_t::gcd_dc() {
//     int g = -1;
//     for (int j=0; j<vd.size(); j++) {
//         if (vd[j] != 0) {
//             g = (g == -1) ? abs(vd[j]) : gcd(g, abs(vd[j]));
//             if (g == 1)
//                 return 1;
//         }
//     }
//     for (int j=0; j<vc.size(); j++) {
//         if (vc[j] != 0) {
//             g = (g == -1) ? abs(vc[j]) : gcd(g, abs(vc[j]));
//             if (g == 1)
//                 return 1;
//         }
//     }
//     return g;
// }

// void row_t::canonicalize() {
//     int div = gcd_dc();
//     if (div > 1) {
//         for (int j=0; j<vd.size(); j++)
//             vd[j] /= div;
//         for (int j=0; j<vc.size(); j++)
//             vc[j] /= div;
//     }
// }

// /////////////////////////////////////////////////////////////////////////////////////////

// std::vector<std::vector<int>>
// fourier_motzkin_kernel(const std::vector<std::vector<int>>& A) 
// {
//     const size_t n = A.size(), m = A.front().size();
//     // cout << "n:"<<n<<" m:"<<m<<endl;
//     std::list<row_t> DC;
//     for (size_t i=0; i<A.size(); i++) {
//         row_t r(n, m);
//         r.vd[i] = 1;
//         r.add_vec(A[i]);
//         DC.emplace_back(std::move(r));
//     }
//     // Find the null space
//     for (int j=0; j<m; j++) {
//         // cout << "j:"<<j<<" DC:"<<DC.size()<<endl;
//         std::list<row_t> NZ;
//         for (auto iter = DC.begin(); iter != DC.end(); ) {
//             int a_ij = iter->vc[j];
//             if (a_ij != 0) {
//                 NZ.splice(NZ.begin(), DC, iter++);
//             }
//             else ++iter;
//         }
//         if (NZ.size() <= 1)
//             continue;

//         auto iter = NZ.begin();
//         ++iter;
//         for (; iter != NZ.end(); ++iter) {
//             int c1 = NZ.front().vc[j];
//             int c2 = iter->vc[j];
//             int k1 = abs(c1), k2 = abs(c2);
//             if (c1 * c2 >= 0)
//                 k2 *= -1;

//             row_t r(k1, *iter, k2, NZ.front());
//             r.canonicalize();
//             DC.emplace_back(std::move(r));
//         }
//     }

//     std::vector<std::vector<int>> out;
//     for (auto&& row : DC) {
//         out.emplace_back(row.vd);
//     }
//     return out;
// }

/////////////////////////////////////////////////////////////////////////////////////////

bool
read_mat(const char* fname,
         std::vector<std::vector<int>>& A) 
{
    A.clear();
    ifstream ifs(fname);
    if (!ifs) return false;
    size_t n, m;
    ifs >> n >> m;
    if (!ifs) return false;

    for (size_t j=0; j<n; j++) {
        std::vector<int> row(m);
        for (size_t i=0; i<m; i++) {
            ifs >> row[i];
            if (!ifs) return false;
        }

        A.emplace_back(std::move(row));
    }
    return true;
}


/////////////////////////////////////////////////////////////////////////////////////////

bool
read_pin(const char* fname, const size_t m,
         std::vector<std::vector<int>>& A) 
{
    A.clear();
    ifstream ifs(fname);
    if (!ifs) return false;
    size_t n;
    ifs >> n;
    if (!ifs) return false;

    for (size_t j=0; j<n; j++) {
        std::vector<int> row(m);
        size_t nnz, index;
        ifs >> nnz;
        if (!ifs) return false;
        
        for (size_t k=0; k<nnz; k++) {
            int value;
            ifs >> value >> index;
            index--; // indices are 1-based in PIN format
            if (!ifs || index>=m) return false;

            row[index] = value; 
        }
        A.emplace_back(std::move(row));
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int>>
transpose(const std::vector<std::vector<int>>& A)
{
    const size_t n = A.size(), m = A.front().size();
    std::vector<std::vector<int>> AT(m);
    for (size_t j=0; j<m; j++)
        AT[j].resize(n);

    for (size_t i=0; i<n; i++)
        for (size_t j=0; j<m; j++)
            AT[j][i] = A[i][j];

    return AT;
}

/////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int>>
make_symmetric(const std::vector<std::vector<int>>& A)
{
    const size_t n = A.size(), m = A.front().size();
    std::vector<std::vector<int>> S(2*n);

    for (size_t i=0; i<n; i++) {
        S[2*i  ].resize(m);
        S[2*i+1].resize(m);
        for (size_t j=0; j<m; j++) {
            S[2*i  ][j] =  A[i][j];
            S[2*i+1][j] = -A[i][j];
        }
    }
    return S;
}

/////////////////////////////////////////////////////////////////////////////////////////

// for each column j<m find the first row i s.t. A[0..i-1, j] <= 0, and A[i,j] > 0
std::vector<size_t>
step_for_negative_removal(const std::vector<std::vector<int>>& A) {
    const size_t n = A.size(), m = A.front().size();
    std::vector<size_t> red(m, 0);

    for (size_t j=0; j<m; j++) {
        for (size_t i=0; i<n; i++) {
            if (A[i][j] != 0)//> 0)
                break;
            else
                red[j]++;
        }
    }

    return red;
}

/////////////////////////////////////////////////////////////////////////////////////////

// replace row i0 with a linear combination of both row i0 and ik s.t. A[i0,j] becomes 0
// Perform integral elementary row operations on row i0 using row ik
// such that the entry A[i0, j] is zeroed
void annul_row(std::vector<std::vector<int>>& A, 
               size_t i0, size_t ik, size_t j) 
{
    const size_t m = A.front().size();
    assert(A[i0][j] != 0 && A[ik][j] != 0);
    int mult_ik = abs(A[i0][j]);
    int mult_i0 = abs(A[ik][j]);
    int gcd_i0k = gcd(mult_i0, mult_ik);
    mult_ik /= gcd_i0k;
    mult_i0 /= gcd_i0k;
    // ensure opposite signs
    if (A[i0][j] * A[ik][j] > 0)
        mult_ik = -mult_ik;
    // Sum row ik to row i0, and verify that B[i0,j] is 0
    for (size_t c=0; c<m; c++) 
        A[i0][c] = mult_ik * A[ik][c] + mult_i0 * A[i0][c];
    assert(A[i0][j] == 0);
    reduce_gcd(A[i0]);
}

// perform integral elementary column operations on column j0 using column jk
// such that the entry A[i, j0] is zeroed. Perform the same operations on matrix E
void annul_columns(std::vector<std::vector<int>>& A, 
                   std::vector<std::vector<int>>& E,
                   size_t j0, size_t jk, size_t i) 
{
    assert(A[i][j0] != 0 && A[i][jk] != 0);
    int mult_jk = abs(A[i][j0]);
    int mult_j0 = abs(A[i][jk]);
    int gcd_j0k = gcd(mult_jk, mult_j0);
    mult_jk /= gcd_j0k;
    mult_j0 /= gcd_j0k;
    // ensure opposing signs
    if (A[i][j0] * A[i][jk] > 0)
        mult_jk = -mult_jk;
    // Sum columns and verify that A[i][j0] is 0
    for (size_t r=0; r<A.size(); r++) 
        A[r][j0] = mult_jk * A[r][jk] + mult_j0 * A[r][j0];
    for (size_t r=0; r<E.size(); r++)
        E[r][j0] = mult_jk * E[r][jk] + mult_j0 * E[r][j0];
    assert(A[i][j0] == 0);
}

/////////////////////////////////////////////////////////////////////////////////////////

void add_mult_row2(std::vector<int>& dst, int mult1, 
                   const std::vector<int>& row2, int mult2)
{
    assert(dst.size() == row2.size());
    for (size_t j=0; j<dst.size(); j++) {
        // dst[j] = mult1 * dst[j] + mult2 * row2[j];
        dst[j] = add_exact(multiply_exact(mult1, dst[j]),
                           multiply_exact(mult2, row2[j]));
    }
}

// replace row i0 with a linear combination of both row i0 and ik s.t. A[i0,j] becomes 0
// Perform integral elementary row operations on row i0 using row ik
// such that the entry A[i0, j] is zeroed
void annul_row_HU(std::vector<std::vector<int>>& A,
                  std::vector<std::vector<int>>& U, // repeat the same operations on U
                  size_t i0, size_t ik, size_t j) 
{
    const size_t m = A.front().size();
    assert(A[i0][j] != 0 && A[ik][j] != 0);
    int mult_ik = abs(A[i0][j]);
    int mult_i0 = abs(A[ik][j]);
    int gcd_i0k = gcd(mult_i0, mult_ik);
    mult_ik /= gcd_i0k;
    mult_i0 /= gcd_i0k;
    // ensure opposite signs
    if (sign3(A[i0][j]) * sign3(A[ik][j]) > 0)
        mult_ik = -mult_ik;
    // Sum row ik to row i0, and verify that B[i0,j] is 0
    cout << "  * annul_row_HU: "<<mult_i0<<"*["<<i0<<"] + "<<mult_ik<<"*["<<ik<<"]" << endl;
    add_mult_row2(A[i0], mult_i0, A[ik], mult_ik);
    add_mult_row2(U[i0], mult_i0, U[ik], mult_ik);
    assert(A[i0][j] == 0);
    // reduce_gcd(A[i0]);
}

// // replace row i0 with a linear combination of both row i0 and ik s.t. A[i0,j] becomes 0
// // Perform integral elementary row operations on row i0 using row ik
// // such that the entry A[i0, j] is zeroed
// void annul_row_H(std::vector<std::vector<int>>& A, 
//                  size_t i0, size_t ik, size_t j) 
// {
//     const size_t m = A.front().size();
//     assert(A[i0][j] != 0 && A[ik][j] != 0);
//     int mult_ik = abs(A[i0][j]);
//     int mult_i0 = abs(A[ik][j]);
//     int gcd_i0k = gcd(mult_i0, mult_ik);
//     mult_ik /= gcd_i0k;
//     mult_i0 /= gcd_i0k;
//     // ensure opposite signs
//     if (A[i0][j] * A[ik][j] > 0)
//         mult_ik = -mult_ik;
//     // Sum row ik to row i0, and verify that B[i0,j] is 0
//     for (size_t c=0; c<m; c++) 
//         A[i0][c] = mult_ik * A[ik][c] + mult_i0 * A[i0][c];
//     assert(A[i0][j] == 0);
//     // canonicalize(A[i0]);
// }

// // perform integral elementary column operations on column j0 using column jk
// // such that the entry A[i, j0] is zeroed. 
// void annul_column_H(std::vector<std::vector<int>>& A, 
//                     size_t j0, size_t jk, size_t i) 
// {
//     assert(A[i][j0] != 0 && A[i][jk] != 0);
//     int mult_jk = abs(A[i][j0]);
//     int mult_j0 = abs(A[i][jk]);
//     int gcd_j0k = gcd(mult_jk, mult_j0);
//     mult_jk /= gcd_j0k;
//     mult_j0 /= gcd_j0k;
//     // ensure opposing signs
//     if (A[i][j0] * A[i][jk] > 0)
//         mult_jk = -mult_jk;
//     // Sum columns and verify that A[i][j0] is 0
//     for (size_t r=0; r<A.size(); r++) 
//         A[r][j0] = mult_jk * A[r][jk] + mult_j0 * A[r][j0];
//     assert(A[i][j0] == 0);
// }

/////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int>> 
mat_mult(const std::vector<std::vector<int>>& A, 
         const std::vector<std::vector<int>>& B) 
{
    const size_t nA = A.size(), mA = A[0].size();
    const size_t nB = B.size(), mB = B[0].size();
    assert(mA == nB);

    // initialize with zeroes
    std::vector<std::vector<int>> C(nA);
    for (size_t i = 0; i < nA; ++i)
        C[i].resize(mB);

    // Multiply
    for (size_t i = 0; i < nA; ++i) {
        for (size_t j = 0; j < mB; ++j) {
            for (size_t k = 0; k < nB; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C; // C = A * B
}

/////////////////////////////////////////////////////////////////////////////////////////

void row_echelon_form(std::vector<std::vector<int>>& A, 
                      std::vector<size_t> *leading_cols) 
{
    if (leading_cols) leading_cols->clear();
    const size_t n = A.size();//, m = A.front().size();
    for (size_t k=0; k<n; k++) {
        int i_max = k; // Find the k-th pivot row
        for (size_t i=k+1; i<n; i++) {
            if (count_nonzero(A[i_max]) == 0 || 
                (count_nonzero(A[i]) > 0 && leading(A[i]) < leading(A[i_max])))
                i_max = i;
        }
        //  Move the pivot row in position k
        if (k != i_max)
            std::swap(A[k], A[i_max]);
        reduce_gcd(A[k]);
        // Annull column j0 to all the rows above and below the pivot (row k)
        size_t j0 = leading(A[k]);
        if (j0 != size_t(-1)) {
            if (leading_cols) leading_cols->push_back(j0);
            for (size_t i=0; i<n; i++) {
                if (i!=k && A[i][j0] != 0)
                    annul_row(A, i, k, j0);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

void row_footprint_form(std::vector<std::vector<int>>& A, 
                        std::vector<size_t> *leading_cols,
                        std::vector<size_t> *trailing_cols)
{
    if (leading_cols) leading_cols->clear();
    if (trailing_cols) trailing_cols->clear();
    const size_t n = A.size();//, m = A.front().size();
    for (size_t k=0; k<n; k++) {
        int i_max = k; // Find the k-th pivot row
        for (size_t i=k+1; i<n; i++) {
            if (count_nonzero(A[i_max]) == 0 || 
                (count_nonzero(A[i]) > 0 && leading(A[i]) < leading(A[i_max])))
                i_max = i;
        }
        //  Move the pivot row in position k
        if (k != i_max)
            std::swap(A[k], A[i_max]);
        reduce_gcd(A[k]);
        // Annull column j0 to all the rows above and below the pivot (row k)
        size_t j0 = leading(A[k]);
        if (j0 != size_t(-1)) {
            if (leading_cols) leading_cols->push_back(j0);
            for (size_t i=k+1; i<n; i++) {
                if (A[i][j0] != 0)
                    annul_row(A, i, k, j0);
            }
        }
    }
    // Step 2: Find row-trailing entries and annul all entries above each of them. 
    for (ssize_t k=n-1; k>=0; k--) {
        // Annul the last column of B[k] to all the rows above the pivot row k
        size_t jN = trailing(A[k]);
        if (jN != size_t(-1)) {
            if (trailing_cols) trailing_cols->push_back(jN);
            for (ssize_t i=k-1; i>=0; i--)
                if (A[i][jN] != 0)
                    annul_row(A, i, k, jN);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

// irank is the sum of the row spans in rrff form
size_t irank(const std::vector<std::vector<int>>& A) {
    auto rrffA = A;
    row_footprint_form(rrffA);

    size_t v = 0;
    for (const auto& row : rrffA) {
        // cout << "  "<<trailing(row)<<" - "<<leading(row)<<endl;
        v += trailing(row) - leading(row) + 1; // row span
    }
    return v;
}

/////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int>>
nullspace_from_hnf(const std::vector<std::vector<int>>& A,
                   const std::vector<size_t>& pivot_cols) 
{
    const size_t m = A.front().size();
    // lcm of all pivot_cols
    int lcm_val = 1;
    for (size_t i=0; i<pivot_cols.size(); i++)
        lcm_val = lcm(lcm_val, A[i][pivot_cols[i]]);

    std::vector<std::vector<int>> ns_basis;
    // for all free variables
    for (size_t j=0; j<m; j++) {
        if (std::find(pivot_cols.begin(), pivot_cols.end(), j) != pivot_cols.end())
            continue; // j is not a free variable

        std::vector<int> row(m);
        row[j] = lcm_val;

        for (size_t i=0; i<pivot_cols.size(); i++)
            row[pivot_cols[i]] -= A[i][j] * (lcm_val / A[i][pivot_cols[i]]);

        reduce_gcd(row);        
        ns_basis.emplace_back(std::move(row));
    }
    return ns_basis;
}

/////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int>>
nullspace(std::vector<std::vector<int>>& A) {
    std::vector<size_t> pivot_cols;
    row_echelon_form(A, &pivot_cols); // bring A in REF

    return nullspace_from_hnf(A, pivot_cols);
}

/////////////////////////////////////////////////////////////////////////////////////////

void swap_columns(std::vector<std::vector<int>>& A, size_t j1, size_t j2) {
    assert(j1 < A.front().size() && j2 < A.front().size());
    for (size_t r=0; r<A.size(); r++)
        std::swap(A[r][j1], A[r][j2]);
}

void negate_row(std::vector<int>& A) {
    for (size_t j=0; j<A.size(); j++)
        A[j] = -A[j];
}

void add_mult_row(std::vector<int>& dst, int mult, 
                  const std::vector<int>& addend)
{
    for (size_t j=0; j<dst.size(); j++)
        dst[j] = add_exact(dst[j], multiply_exact(mult, addend[j]));
}

void add_mult_column(std::vector<std::vector<int>>& A,
                     const size_t j0, const int mult, const size_t jk) 
{
    for (size_t r=0; r<A.size(); r++)
        A[r][j0] = add_exact(A[r][j0], multiply_exact(mult, A[r][jk]));
}

int column_lcm(const std::vector<std::vector<int>>& A, size_t j) {
    int l = 0;
    for (size_t r=0; r<A.size(); r++)
        if (A[r][j])
            l = (l == 0) ? abs(A[r][j]) : lcm(l, abs(A[r][j]));
    return l;
}

int column_gcd(const std::vector<std::vector<int>>& A, size_t j) {
    int g = 0;
    for (size_t r=0; r<A.size(); r++)
        if (A[r][j])
            g = gcd(g, abs(A[r][j]));
    return g;
}

// Reduce gcd of @row. If gcd is not 1, then
// multiply all rows of U by g except row i0
void reduce_gcd_HU(std::vector<int>& row, 
                   std::vector<std::vector<int>>& U, size_t i0) 
{
    int g = reduce_gcd(row);
    if (g > 1) {
        for (size_t i=0; i<U.size(); i++) {
            if (i != i0) {
                for (size_t j=0; j<U[i].size(); j++)
                    U[i][j] = multiply_exact(g, U[i][j]);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////

// rewrite A in block form using row & column swaps to isolate an identity block, as:
//   [ A' 0 ]
//   [ 0  I ]
// where I has size len(col_pivots) * len(col_pivots)
void move_identity_rows_to_right(std::vector<std::vector<int>>& H,
                                 std::vector<size_t>& col_pivots) 
{
    const size_t nr = H.size();
    const size_t nc = H.front().size();
    size_t num_idents = 0;
    col_pivots.clear();
    // identify and separate column identities
    for (ssize_t j=nc-1; j>=0; j--) {
        bool is_identity = true;
        size_t id_row = -1, count_ones = 0;
        for (size_t i=0; i<nr; i++) {
            if (count_ones==0 && abs(H[i][j]) == 1 && i<nr-num_idents) {
                id_row = i;
                count_ones++;
            }
            else if (H[i][j] != 0) {
                // multiple nnz, not a 1 or nnz in the last @num_idents rows
                is_identity = false;
                break;
            }
        }
        if (is_identity && count_ones==1) { // move identity row to bottom-right
            cout << "Identity row/column: "<<id_row<<" moved to j="<<j<<endl;
            std::swap(H[id_row], H[nr-num_idents-1]);
            swap_columns(H, j, nc-num_idents-1);
            col_pivots.push_back(j);
            num_idents += 1;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int>>
identity(size_t n, size_t m) {
    std::vector<std::vector<int>> I(n);
    for (size_t i=0; i<n; i++) {
        I[i].resize(m);
        for (int j=0; j<m; j++)
            I[i][j] = (i==j) ? 1 : 0;
    }
    return I;
}

/////////////////////////////////////////////////////////////////////////////////////////

template<class iter>
int gcd_seq(iter it, const iter end) {
    int g = 0;
    while (it != end) {
        g = gcd(g, abs(*it));
        ++it;
    }
    return g;
}

/////////////////////////////////////////////////////////////////////////////////////////

void print_mat_HU(const std::vector<std::vector<int>>& H,
                  const std::vector<std::vector<int>>& U,
                  size_t iden, size_t i0)
{
    const size_t nc = H.front().size();
    const size_t nr = H.size();
    const size_t ncU = U.front().size();
    const size_t nrU = U.size();
    const size_t m = min(nc, nr) - iden;
    size_t iid, bri;

    // print H   ┼─│║═╬
    cout << "H:\n";
    iid = nc-iden-1;
    bri = nrU-m+1; // bottom-right identity
    for (size_t i=0; i<nr; i++) {
        for (size_t j=0; j<nc; j++) {
            cout << setw(3) << H[i][j];
            if (j==i0 || j==iid || j==bri) cout << "│";
            else cout << " ";
        }
        cout << endl;
        if (i==i0 || i==m-1) {
            for (size_t j=0; j<nc; j++) {
                cout << "───";
                if (j==i0 || j==iid || j==bri) cout << "┼"; else cout << "─";
            }
            cout << endl;
        }
    }
    // print U
    cout << "U:\n";
    for (size_t i=0; i<nrU; i++) {
        for (size_t j=0; j<ncU; j++) {
            cout << setw(3) << U.at(i).at(j) ;
            if (j==bri) cout << "│";
            else cout << " ";
        }
        cout << endl;
        if (i==bri) {
            for (size_t j=0; j<ncU; j++) {
                cout << "───";
                if (j==bri) cout << "┼"; else cout << "─";
            }
            cout << endl;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

// compute an integral kernel of kernel of matrix A using the
// problem rewriting defined by Hemmecke (used by 4ti2)
// Code derived from 4ti2/src/zsolve/Lattice.hpp (GPL-v2.0)
// for benchmark and comparisons
bool
integral_kernel_old(const std::vector<std::vector<int>>& A,
                    std::vector<std::vector<int>>& result,
                    bool verbose)
{
    // cout << "A:\n"; print_mat(A);
    const size_t ncA = A.front().size();
    const size_t nrA = A.size();
    // Canonicalize rows
    std::vector<std::vector<int>> H = A;
    for (size_t i=0; i<nrA; i++)
        reduce_gcd(H[i]);
    // Build the nullspace U by reducing H while keeping track of the operations
    std::vector<size_t> idents_pivots;
    move_identity_rows_to_right(H, idents_pivots);
    if (verbose) {
        cout << "Moving "<<idents_pivots.size()<<" identity rows/cols to the bottom-right.\n";
        cout << "H:\n"; print_mat(H); cout << endl;
    }
    size_t nc = ncA - idents_pivots.size();
    size_t nr = nrA - idents_pivots.size();
    if (nc==0) 
        return false;
    
    const std::vector<std::vector<int>> C = H;
    std::vector<std::vector<int>> U = identity(nc, nc);
    // std::vector<std::vector<int>> P = identity(nrA, nrA);
    const size_t m = min(nc, nr);
    if (verbose) {
        cout << "ncA="<<ncA<<" nrA="<<nrA<<"  nc="<<nc<<" nr="<<nr<<"  m="<<m<<endl;
        cout << "shape(A,H) = "<<nrA<<"*"<<ncA<<"   shape(U) = "<<nc<<"*"<<nc<<endl;
    }
    size_t i = 0;
    while (i < m) {
        if (verbose) { 
            cout << "\n\nStart of new iteration for i="<<i<<endl;
            print_mat_HU(H, U, idents_pivots.size(), i-1); cout << endl;
        }
        // find pivot row with smallest gcd and move it to row index i
        size_t rpivot_i = i;
        int rpivot_val = -1;
        for (int i2=i; i2<nr; i2++) {
            int rval = gcd_seq(H[i2].begin()+i, H[i2].begin()+nc);
            // cout << "row gcd("<<i2<<", "<<i<<":"<<nc<<") = "<<rval<<endl;
            if (rval > 0 && (rpivot_val == -1 || rval < rpivot_val)) {
                rpivot_i = i2;
                rpivot_val = rval;
            }
        }
        if (rpivot_val == -1)
            break; // no more pivots

        if (verbose) cout << "Swap row pivot="<<rpivot_i<<" -> "<<i<<endl;
        std::swap(H[i], H[rpivot_i]);
        // std::swap(U[i], U[rpivot_i]); //TODO: remove

        // cout << "row pivot="<<rpivot_i<<" -> "<<i<<"     rpivot_val="<<rpivot_val<<endl; // TODO: remove
    
        bool reduced = true;
        while (reduced) {
            reduced = false;
            // find a reducer column j with smallest H[i,j]
            size_t cpivot_j = -1;
            for (size_t j2=i; j2<nc; j2++) {
                if (abs(H[i][j2]) > 0 && 
                    (cpivot_j == -1 || abs(H[i][j2]) < abs(H[i][cpivot_j])))
                {
                    cpivot_j = j2;
                }
            }
            assert(cpivot_j != -1);
            // cout << " pivot i:"<<i<<" is column:"<<cpivot_j<<"   rpivot_val="<<abs(H[i][cpivot_j])<<endl;
            if (verbose)
                cout << "Column "<<cpivot_j<<" will reduce row "<<i
                    << " of other columns. Value="<<abs(H[i][cpivot_j])<<endl;

            // reduce other non-pivot columns
            for (size_t j=0; j<nc; j++) {
                if (j != cpivot_j) {
                    if (H[i][j] != 0) {
                        if (verbose)
                            cout << " column "<<cpivot_j<<" annuls column "<<j<<" at row "<<i<<endl;
                        // assert(0 == (abs(H[i][j]) %  abs(H[i][cpivot_j])));
                        // annul_columns(H, U, j, cpivot_j, i);
                        int mult =  -H[i][j] / H[i][cpivot_j];

                        // int v = H[i][j], r=abs(H[i][cpivot_j]), mult=0;
                        // cout <<"H[i][j]="<<H[i][j]<<" H[i][cpivot_j]="<<H[i][cpivot_j]<<endl;
                        // while (v < 0) { v += r; mult++; }
                        // while (v >= r) { v -= r; mult--; }
                        // if (H[i][cpivot_j] < 0) mult = -mult;

                        if (mult != 0) {
                            if (verbose)
                                cout <<"H[i][j]="<<H[i][j]<<" bv="<<abs(H[i][cpivot_j])
                                     <<" mult="<<mult<<endl;
                            // cout << "  + add column "<<cpivot_j<<" to "<<j<<endl;
                            add_mult_column(H, j, mult, cpivot_j);
                            add_mult_column(U, j, mult, cpivot_j);
                            assert(abs(H[i][j]) < abs(H[i][cpivot_j]));
                            reduced = true;
                        }
                    }
                }
            }
            if (verbose) cout << "Swap column pivot:"<<cpivot_j<<" -> "<<i<<endl;
            swap_columns(H, i, cpivot_j);
            swap_columns(U, i, cpivot_j);
        }
        i++;
    }
    if (verbose) {
        cout << "\n\nFinal matrices:"<<endl;
        print_mat_HU(H, U, idents_pivots.size(), m-1); cout << endl;
    }
        
    // build the integral kernel from U columns
    result.clear(); // nc * ncA
    while (i < nc) {
        std::vector<int> vec(ncA);
        for (size_t j=0; j<nc; j++) // non-identity columns
            vec[j] = U[j][i];
        for (size_t jj=0; jj<idents_pivots.size(); jj++) { // identity columns
            int mult = 0;
            for (size_t k=0; k<nc; k++)
                mult -= vec[k] * C[nr+jj][k];
            vec[nc+jj] = mult * H[nr+jj][nc+jj];
        }
        // reduce_gcd(vec); // NO: this breaks the lin.eq.system
        result.emplace_back(std::move(vec));
        i++;
    }
    // reorder identities
    for (ssize_t c=idents_pivots.size()-1; c>=0; c--)
        swap_columns(result, idents_pivots[c], ncA-c-1);

    // if (verbose) {
    //     std::vector<std::vector<int>> U2 = U; // nc * nc
    //     // U2.resize(m); // m * nc
    //     // U2 = transpose(U2); // nc * m
    //     std::vector<std::vector<int>> A2 = C; // nrA * ncA
    //     cout << "A2:\n"; print_mat(A2); cout << endl;
    //     // A2.resize(m); // m * ncA
    //     for (auto&& row : A2) row.resize(nc); // m * nc
    //     A2 = transpose(A2); // nc * m
    //     cout << "U: "<<U2.size()<<" * "<<U2.front().size()<<endl;
    //     cout << "A: "<<A2.size()<<" * "<<A2.front().size()<<endl;
    //     cout << "H: "<<H.size()<<" * "<<H.front().size()<<endl;
    //     std::vector<std::vector<int>> UA = mat_mult(U2, A2);
    //     UA = transpose(UA);
    //     cout << "U * A:\n"; print_mat(UA); cout << endl;
    // }

    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Hermite normal form
/////////////////////////////////////////////////////////////////////////////////////////

static void
print_mat_hnf(const std::vector<std::vector<int>>& H, // n*m
              const std::vector<std::vector<int>>& U, // n*n
              const std::vector<size_t>& leading_cols,
              const size_t j0) 
{
    const size_t n = H.size();
    const size_t m = H.front().size();

    auto HH = H;
    auto UU = U;
    std::vector<size_t> indices(m);
    // for (size_t j=0; j<m; j++)
    //     indices[j] = j;

    // from H to HH columns
    for (size_t k=0, j=0; k<m; k++) {
        if(k < leading_cols.size()) {
            indices[k] = leading_cols[k];
        }
        else {
            // find next j not in leading_cols[]
            bool found = false;
            while (!found) {
                if (std::find(leading_cols.begin(), leading_cols.end(), j) == leading_cols.end()) {
                    found = true;
                }
                else j++;
            }
            indices[k] = j++;
        }
    }
    for (size_t i=0; i<n; i++) {
        for (size_t j=0; j<m; j++)
            HH[i][j] = H[i][indices[j]];
        for (size_t j=0; j<n; j++)
            UU[i][j] = U[i][j<m ? indices[j] : j];
    }

    // for (size_t k=0; k<leading_cols.size(); k++) {
    // // for (ssize_t k=leading_cols.size()-1; k>=0; k--) {
    //     // cout << "swap columns: "<<indices[k]<<" <> "<<leading_cols[k]<<endl;
    //     swap_columns(HH, indices[k], leading_cols[k]);
    //     // swap_columns(UU, indices[k], leading_cols[k]); // columns are over m, UU is n*n
    //     std::swap(indices[k], indices[leading_cols[k]]);
    // }
    cout << "leading_cols: ";
    for (size_t k=0; k<leading_cols.size(); k++)
        cout << leading_cols[k] << " ";
    cout << endl;

    cout << "indices:      ";
    for (size_t i : indices) cout << i << " ";
    cout << endl;

    // print H   ┼─│║═╬
    cout << "H:\n";
    for (size_t i=0; i<n; i++) {
        for (size_t j=0; j<m; j++) {
            cout << setw(3) << HH[i][j];
            if (j==j0) cout << "│";
            else cout << " ";
        }
        cout << endl;
        if (i==j0 || i==m-1) {
            for (size_t j=0; j<m; j++) {
                cout << "───";
                if (j==j0) cout << "┼"; else cout << "─";
            }
            cout << endl;
        }
    }
    // print U
    cout << "U:\n";
    for (size_t i=0; i<n; i++) {
        for (size_t j=0; j<n; j++) {
            cout << setw(3) << UU.at(i).at(j);
            if (j==j0) cout << "│";
            else cout << " ";
        }
        cout << endl;
        if (i==j0) {
            for (size_t j=0; j<n; j++) {
                cout << "───";
                if (j==j0) cout << "┼"; else cout << "─";
            }
            cout << endl;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

static std::vector<int>
hnf_preanalysis(const std::vector<std::vector<int>>& H, size_t k) 
{
    const size_t n = H.size(), m = H.front().size();
    std::vector<int> col_weights(m);

    for (size_t i=k; i<n; i++) {
        size_t nnz = 0, nnzid = size_t(-1);
        for (size_t j=k; j<m; j++) {
            if (H[i][j] != 0) {
                ++nnz;
                nnzid = j;
            }
        }
        if (nnz == 1) {
            col_weights[nnzid]++;
            // cout << "Column "<<nnzid<<" has value "<<col_weights[nnzid]<<"."<<endl;
        }
    }

    return col_weights;
}

/////////////////////////////////////////////////////////////////////////////////////////

static inline bool
hnf_next_column_heur(const std::vector<std::vector<int>>& H, 
                     const std::vector<bool>& pivot_cols,
                     const size_t k, size_t &j,
                     const bool forced_pivot_order,
                     size_t &fpo)
{
    const size_t n = H.size(), m = H.front().size();

    if (forced_pivot_order) {
        // we already have a pivot order. Just search the next non-zero column
        // and store the number in fpo (going in decreasing order)
        while (fpo>0) // forced pivot order index
        {
            fpo--;
            assert(fpo < m);
            assert(!pivot_cols[fpo]);
            bool col_nnz = true;
            for (size_t i=k; i<n && col_nnz; i++) {
                if (H[i][fpo] != 0) {
                    col_nnz = false;
                }
            }
            if (!col_nnz) {
                j = fpo;
                return true;
            }
        }
        return false;
    }

    std::vector<int> cw = hnf_preanalysis(H, k);

    j = -1;
    std::tuple<int, int, int> cpivot_val(0,0,0);
    for (size_t j2=0; j2<m; j2++) {
        // find a pivot column j2 with smallest gcd.
        if (!pivot_cols[j2]) { 

            std::tuple<int, int, int> w(0,0,0);
            for (size_t i=k; i<n; i++) {
                if (H[i][j2] != 0) {
                    std::get<0>(w) = (cw[j2] == 1);
                    std::get<1>(w) = gcd(std::get<1>(w), abs(H[i][j2]));
                    std::get<2>(w) += 1;
                    // std::get<2>(w) += abs(H[i][j2]);
                }
            }

            // cout << "  j2:"<<j2<<"/"<<m<<"  w:"<<std::get<0>(w)<<","<<std::get<1>(w)<<","<<std::get<2>(w)<<endl;
            if (std::get<1>(w) > 0) {
                if (j == size_t(-1) || w < cpivot_val)
                {
                    j = j2;
                    cpivot_val = w;
                }
            }
        }
    }
    // cout << "Next column: "<<j<<"   cpivot_val="<<std::get<0>(cpivot_val)<<","
    //      <<std::get<1>(cpivot_val)<<","<<std::get<2>(cpivot_val)<<","<<endl;

    return (j != -1);
}

/////////////////////////////////////////////////////////////////////////////////////////

static inline void
hnf_next_pivot(const std::vector<std::vector<int>>& H, 
               const std::vector<bool>& pivot_cols,
               const size_t k, const size_t j,
               size_t &i)
{
    const size_t n = H.size(), m = H.front().size();
    i = size_t(-1);

    std::pair<int,int> rpivot_val(0,0);
    for (int i0=k; i0<n; i0++) {

        std::pair<int,int> w(abs(H[i0][j]),0);
        if (w.first > 0) {
             for (size_t j2=0; j2<m; j2++) {
                if (!pivot_cols[j2]) { 
                    if (H[i0][j2] != 0)
                        w.second++;
                }
             }
        }

        if (w.first > 0 && (i == size_t(-1) || w < rpivot_val)) {
            i = i0;
            rpivot_val = w;
        }
    }

    assert(i != size_t(-1));
}


/////////////////////////////////////////////////////////////////////////////////////////

// // get the gcd for column j for the rows k..n
// int column_gcd_range(const std::vector<std::vector<int>>& A, 
//                      size_t j, size_t k) 
// {
//     const size_t n = A.size();

//     int g = 0;
//     for (size_t i=k; i<n; i++) {
//         if (A[i][j] != 0) {
//             g = gcd(g, abs(A[i][j]));
//         }
//     }
//     return g; 
// }

/////////////////////////////////////////////////////////////////////////////////////////

// Compute row-style Hermite normal form ( U * A = H ) where:
//  - H is upper triangular (after swapping the leading_cols[])
//  - U is unimodular (invertible over the integer domain)
//  - leading coefficients of non-zero rows are positive
//  - elements below pivots are zero
//  - elements above pivots are non-negative and strictly smaller than the pivot
//  - pivot column for non-zero row i is stored in leading_cols[i]
// See Algorithm 2.4.5, "A course in computational algebraic number theory", Cohen H.
void 
hermite_normal_form(const std::vector<std::vector<int>>& A, 
                    std::vector<std::vector<int>>& H,
                    std::vector<std::vector<int>>& U,
                    std::vector<size_t> *leading_cols,
                    const bool forced_pivot_order,
                    bool verbose) 
{
    std::vector<size_t> _lc;
    if (!leading_cols) leading_cols = &_lc;
    leading_cols->clear();
    const size_t n = A.size();
    const size_t m = A.front().size();
    H = A; 
    U = identity(n, n);
    if (n==0)
        return;

    std::vector<bool> pivot_cols(m);

    size_t k, fpo=m;
    for (k=0; k<n; k++) {
        if (verbose) {
            cout << "\n\nStart of new iteration for k="<<k<<endl;
            print_mat_hnf(H, U, *leading_cols, k-1); cout << endl;
        }
        // Find the pivot column
        size_t j;
        if (!hnf_next_column_heur(H, pivot_cols, k, j, forced_pivot_order, fpo))
            break; // no more pivots

        if (verbose) {
            cout << "Pivot column="<<j<<endl;
        }

        leading_cols->push_back(j);
        pivot_cols[j] = true;

        // Euclidean step.
        // Try to reduce rows as much as possibile, only by performing sums
        // with the pivot rows. Never multiply the pivot row.
        bool can_improve_pivot;
        do {
            can_improve_pivot = false;

            // find pivot row with smallest value on column j
            size_t rpivot_k;
            hnf_next_pivot(H, pivot_cols, k, j, rpivot_k);
            if (verbose) {
                cout << "Pivot row="<<rpivot_k<<" col="<<j<<"   H[i,j]="<<H[rpivot_k][j]<<endl;
            }

            // Move the pivot row in position k
            std::swap(H[k], H[rpivot_k]);
            std::swap(U[k], U[rpivot_k]); // TODO: check

            assert(j != size_t(-1) && H[k][j] != 0);
            // make leading positive
            if (H[k][j] < 0) {
                if (verbose) cout << "Making pivot row' leading value positive." << endl;
                negate_row(H[k]);
                negate_row(U[k]);
            }

            // reduce pivot j of all other rows 0..k-1, k+1..n
            for (size_t i0=0; i0<n; i0++) {
                if (i0 != k) {
                    // TODO: Test which combination is better.
                    int mult = -H[i0][j] / H[k][j];
                    // int v = H[i0][j], r=abs(H[k][j]), mult=0;
                    // while (v < 0) { v += r; mult++; }
                    // while (v >= r) { v -= r; mult--; }
                    // cout <<"H[i0][j]="<<H[i0][j]<<" H[k][j]="<<H[k][j]<<" mult="<<mult<<" v="<<v<<endl;
                    
                    if (mult != 0) {
                        // cout << "  + add row "<<k<<" to "<<i0<<"  mult="<<mult<<endl;
                        // subtract row k to non-pivot row i0
                        add_mult_row(H[i0], mult, H[k]); // H[i0] += mult * H[k]
                        add_mult_row(U[i0], mult, U[k]);
                        // reduce_gcd_HU(H[i0], U, i0);
                    }
                    if (!(0 <= abs(H[i0][j]) && abs(H[i0][j]) < H[k][j])) {
                        cout << "i0:"<<i0<<" k:"<<k<<" j:"<<j<<" mult:"<<mult<<endl;
                        cout << "H:" << endl; print_mat(H); cout << flush << endl;
                        throw std::exception();
                    }
                    assert(0 <= abs(H[i0][j]) && abs(H[i0][j]) < H[k][j]);

                    // Check if H[i0][j] is a better pivot than H[k][j] after subtraction
                    // (To be a better pivot it needs to be below pivot row k)
                    if (0 != H[i0][j] && abs(H[i0][j]) < H[k][j] && i0 > k) {
                        can_improve_pivot = true;
                    }
                }
            }

            if (can_improve_pivot && verbose) {
                cout << "\n\nNew row-pivot round for k="<<k<<endl;
                print_mat_hnf(H, U, *leading_cols, k); cout << endl;
            }                
        }
        while (can_improve_pivot);
    }

    // cout << "\n\nFinal matrices for Hermite Normal Form:"<<endl;
    // print_mat_HU(transpose(H), transpose(U), 0, m); cout << endl;
    // std::vector<std::vector<int>> H2 = H;
    // for (size_t k=0; k<leading_cols->size(); k++)
    //     swap_columns(H2, (*leading_cols)[k], k);
    // cout << "Hermite normal form:\n";
    // print_mat(H2); cout << endl;

// #ifndef NDEBUG
//     std::vector<std::vector<int>> H2 = mat_mult(U, A);
//     for (size_t i=0; i<n; i++) 
//         reduce_gcd(H2[i]);
//     // if (H != H2) {
//         cout << "A:" << endl; print_mat(A); cout << endl;
//         cout << "U:" << endl; print_mat(U); cout << endl;
//         cout << "H:" << endl; print_mat(H); cout << endl;
//         cout << "H2:" << endl; print_mat(H2); cout << endl;
//         // assert(false);
//     // }
// #endif

    return;
}

/////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int>>
integral_kernel_Zgens(const std::vector<std::vector<int>>& A,  
                      std::vector<size_t>& leading_cols, bool verbose)
{
    std::vector<std::vector<int>> H, U, UA, basisZkerA;
    leading_cols.clear();
    const size_t n = A.size();

    auto trA = transpose(A);
    // reduce gcd's of A^T columns
    for (size_t j=0; j<n; j++)
        reduce_column_gcd(trA, j);

    // compute: H = U * transpose(A) * T
    hermite_normal_form(trA, H, U, &leading_cols, false, verbose);
    const size_t rankA = leading_cols.size();

    // Take the last n-rank(A) rows of U as the Z-basis for the integral kernel of A
    for (size_t i=rankA; i<U.size(); i++) {
        basisZkerA.emplace_back(std::move(U[i]));
        reduce_gcd(basisZkerA.back());
    }
    // cout << "basisZkerA(A):" << endl; print_mat(basisZkerA); cout << endl;

    return basisZkerA;
}

/////////////////////////////////////////////////////////////////////////////////////////

void hnf_scores(const std::vector<std::vector<int>>& H) {
    const size_t n = H.size();
    if (n==0)
        return;
    const size_t m = H.front().size();
    size_t S1 = 0, S2 = 0, S3 = 0, S4 = 0;

    for (size_t j=0; j<m; j++) {
        size_t n_pos=0, c_pos=0, n_neg=0, c_neg=0;
        for (size_t i=0; i<n; i++) {
            n_pos += H[i][j] > 0;
            c_pos += H[i][j] > 0 ? H[i][j] : 0;
            n_neg += H[i][j] < 0;
            c_neg += H[i][j] < 0 ? -H[i][j] : 0;
        }
        S1 += n_pos + n_neg;
        S2 += (n_pos * n_neg) + n_pos + n_neg;
        S3 += c_pos + c_neg;
        S4 += (c_pos * c_neg) + c_pos + c_neg;
    }

    cout << "S1: " << S1 << endl;
    cout << "S2: " << S2 << endl;
    cout << "S3: " << S3 << endl;
    cout << "S4: " << S4 << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////


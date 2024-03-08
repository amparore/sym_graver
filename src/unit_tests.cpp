#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>

#include "sym_pottier.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////

void unit_test_vcanon() {
    const std::vector<std::vector<int>> mA = {
        { 0, 3, 1, 0, 0},
        { 2, 0, 0, 1, 2},
        { 1,-1,-3,-1,-1},
        { 0, 2, 0, 2, 0},
        { 0, 2, 4, 2, 0},
        { 0, 2, 3, 2, 0},
    };
    const std::vector<std::vector<int>> mB = {
        { 3,-2,-2, 1, 4},
        {-5, 1, 1, 2,-3},
        { 0, 1, 3, 0, 0},
    };
    const size_t num_levels = mA[0].size();
    variable_order vorder(num_levels), pivot_order(num_levels);

    meddly_context ctx(num_levels, vorder, pivot_order);
    size_t meddly_cache_size = 1000000;
    ctx.initialize(meddly_cache_size);
    MEDDLY::dd_edge A = mdd_from_vectors(mA, ctx.forestMDD, false);
    MEDDLY::dd_edge B = mdd_from_vectors(mB, ctx.forestMDD, false);

    cout << "A:\n" << print_mdd(A, vorder) << endl;

    MEDDLY::dd_edge U(ctx.forestMDD), D(ctx.forestMDD);
    VDIVIDE->computeDDEdge(A, 2, U, D);
    cout << "U:\n" << print_mdd(U, vorder) << endl;
    cout << "D:\n" << print_mdd(D, vorder) << endl;

    exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

void unit_test_svectors() {
    const std::vector<std::vector<int>> mA = {
        {  0,  5,  9,  13,  1, 1 }, 
        {  0, -5, -9, -14, -2, 1 }, 
        { -5,  0, -6, -12, -4, 1 }, 
     };
    const std::vector<std::vector<int>> mB = {      
    };
    const size_t num_levels = mA[0].size();
    variable_order vorder(num_levels), pivot_order(num_levels);

    meddly_context ctx(num_levels, vorder, pivot_order);
    size_t meddly_cache_size = 1000000;
    ctx.initialize(meddly_cache_size);
    MEDDLY::dd_edge A = mdd_from_vectors(mA, ctx.forestMDD, false);
    MEDDLY::dd_edge B = mdd_from_vectors(mB, ctx.forestMDD, false);
    B = A;

    cout << "A:\n" << print_mdd(A, vorder) << endl;
    cout << "B:\n" << print_mdd(B, vorder) << endl;
    cout << "----------------------------------\n";

    MEDDLY::dd_edge SVp(ctx.forestMDD);
    MEDDLY::dd_edge SVn(ctx.forestMDD);
    S_VECTORS->computeDDEdge(A, B, true, ab_sum_t::A_PLUS_B, SVS_POS, 6, SVp);
    cout << "A + B:\n" << print_mdd(SVp, vorder) << endl;

    S_VECTORS->computeDDEdge(A, B, true, ab_sum_t::A_MINUS_B, SVS_UNDECIDED, 6, SVn);
    cout << "A - B:\n" << print_mdd(SVn, vorder) << endl;

    exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

void unit_test_reduce() {
    const std::vector<std::vector<int>> mA = {
        { 0, 3, 1 },
        { 0,-3, 1 },
        { 7, 0, 0 },
        { -7,0, 1 },
        // { -2, -2, 2 },
        // { 0, 3, 1, 1 }
     };
    const std::vector<std::vector<int>> mB = {      
        { 0, 4, 1 },
        { 0, 1, 1 },
        { 7, 0, 0 },
        { 0, -1, 0 },
        // { 0, 2, 1, 1 },
        // { 0, 1, 2, 1 },
        // { 2, 1, 1, 0 },
        // { 1, 1, 0, 1 },
        // { 0, 0, 1, 2 },
    };
    const size_t num_levels = mA[0].size();
    variable_order vorder(num_levels), pivot_order(num_levels);

    meddly_context ctx(num_levels, vorder, pivot_order);
    size_t meddly_cache_size = 1000000;
    ctx.initialize(meddly_cache_size);
    MEDDLY::dd_edge A = mdd_from_vectors(mA, ctx.forestMDD, false);
    MEDDLY::dd_edge B = mdd_from_vectors(mB, ctx.forestMDD, false);

    cout << "A:\n" << print_mdd(A, vorder) << endl;
    cout << "B:\n" << print_mdd(B, vorder) << endl;
    cout << "----------------------------------\n";

    MEDDLY::dd_edge C(ctx.forestMDD), S(ctx.forestMDD);
    LEQ_NEQ_SQ_COMPARE->computeDDEdge(A, B, true, true, 0, C);
    LEQ_NEQ_SQ_SUBTRACT->computeDDEdge(A, B, true, true, 0, S);
    cout << "A \\ C:\n" << print_mdd(sym_difference(A, C), vorder) << endl;
    cout << "C:\n" << print_mdd(C, vorder) << endl;
    cout << "S:\n" << print_mdd(S, vorder) << endl;
    cout << "----------------------------------\n";

    MEDDLY::dd_edge I(ctx.forestMDD), R(ctx.forestMDD), D(ctx.forestMDD);
    REDUCE->computeDDEdge(A, B, true, true, sv_sign::SVS_UNDECIDED, cmp_sign::CMP_UNDECIDED, 0, I, D);
    R = sym_difference(A, I);
    cout << "I:\n" << print_mdd(I, vorder) << endl;
    cout << "R:\n" << print_mdd(R, vorder) << endl;
    cout << "D:\n" << print_mdd(D, vorder) << endl;

    // cout << endl;
    // REDUCE3->computeDDEdge(A, B, true, true, sv_sign::SVS_UNDECIDED, cmp_sign::CMP_UNDECIDED, 0, I, R, D);
    // cout << "I:\n" << print_mdd(I, vorder) << endl;
    // cout << "R:\n" << print_mdd(R, vorder) << endl;
    // cout << "D:\n" << print_mdd(D, vorder) << endl;

    exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

void unit_test_sign_canon() {
    const std::vector<std::vector<int>> mA = {
        { 2, 3, 1, 0, 0 },
        {-2,-3,-1, 0, 0 },
        { 5, 4,-1, 0, 0 },
        { 2, 0, 0, 1, 2 },
        { 2, 0,-1, 1,-2 },
    };

    const size_t num_levels = mA[0].size();
    variable_order vorder(num_levels), pivot_order(num_levels);

    meddly_context ctx(num_levels, vorder, pivot_order);
    size_t meddly_cache_size = 1000000;
    ctx.initialize(meddly_cache_size);
    MEDDLY::dd_edge A = mdd_from_vectors(mA, ctx.forestMDD, false);

    cout << "A:\n" << print_mdd(A, vorder) << endl;
    MEDDLY::dd_edge C(ctx.forestMDD);
    SIGN_CANON->computeDDEdge(A, C, false);
    cout << "C:\n" << print_mdd(C, vorder) << endl;

    exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

void unit_test_minimal_supports() {
    const std::vector<std::vector<int>> mA = {
        {  4, 0, 2 },
        {  6, 2, 0 },
        { 12, 4, 0 },
    };

    const size_t num_levels = mA[0].size();
    variable_order vorder(num_levels), pivot_order(num_levels);

    meddly_context ctx(num_levels, vorder, pivot_order);
    size_t meddly_cache_size = 1000000;
    ctx.initialize(meddly_cache_size);
    MEDDLY::dd_edge A = mdd_from_vectors(mA, ctx.forestMDD, false);

    MEDDLY::dd_edge non_minimal_support(A.getForest());
    SUPPORT_INCL_TABLE->get_op(0, true, false)->computeDDEdge(A, A, non_minimal_support, false);
    cout << "A:\n" << print_mdd(A, vorder) << endl;
    cout << "nms:\n" << print_mdd(non_minimal_support, vorder) << endl;
    // A = sym_difference(A, non_minimal_support);
}
    
/////////////////////////////////////////////////////////////////////////////////////////


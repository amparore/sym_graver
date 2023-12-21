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
    VCANON->computeDDEdge(A, 2, U, D);
    cout << "U:\n" << print_mdd(U, vorder) << endl;
    cout << "D:\n" << print_mdd(D, vorder) << endl;

    exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

void unit_test_reduce() {
    const std::vector<std::vector<int>> mA = {
        { 0, 3, 1, 0, 0 },
        { 0,-3, 1, 0, 0 },
        { 2, 0, 0, 1, 2 },
        { 2, 0,-1, 1, 2 },
    };
    const std::vector<std::vector<int>> mB = {
        { 0, 4, 1, 0, 0 },
        { 0, 1, 1, 0, 0 },
        { 2, 0, 0, 1, 2 },
        // { 0, 0, 0, 0, 0 },
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

    MEDDLY::dd_edge U(ctx.forestMDD), R(ctx.forestMDD);
    REDUCE->computeDDEdge(A, B, true, true, 0, U, R);
    cout << "U:\n" << print_mdd(U, vorder) << endl;
    cout << "R:\n" << print_mdd(R, vorder) << endl;

    exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////

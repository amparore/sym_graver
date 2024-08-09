#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>

#include "math_utils.h"
#include "sym_pottier.h"

using namespace std;


/////////////////////////////////////////////////////////////////////////////////////////

void sep(const pottier_params_t& pparams) {
    if (pparams.verbose)
        cout << "----------------------------------------------------------------------" << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge 
sym_union(const MEDDLY::dd_edge& A, const MEDDLY::dd_edge& B) {
    MEDDLY::dd_edge R(A.getForest());
    MEDDLY::apply(MEDDLY::UNION, A, B, R); 
    return R;
}

MEDDLY::dd_edge 
sym_intersection(const MEDDLY::dd_edge& A, const MEDDLY::dd_edge& B) {
    MEDDLY::dd_edge R(A.getForest());
    MEDDLY::apply(MEDDLY::INTERSECTION, A, B, R); 
    return R;
}

MEDDLY::dd_edge 
sym_difference(const MEDDLY::dd_edge& A, const MEDDLY::dd_edge& B) {
    MEDDLY::dd_edge R(A.getForest());
    MEDDLY::apply(MEDDLY::DIFFERENCE, A, B, R); 
    return R;
}

/////////////////////////////////////////////////////////////////////////////////////////

void meddly_context::initialize(size_t meddly_cache_size) 
{
    MEDDLY::initializer_list* L = MEDDLY::defaultInitializerList(0);
    MEDDLY::ct_initializer::setMaxSize(meddly_cache_size);
    MEDDLY::initialize(L);

    // Initialize domain
    MEDDLY::domain *dom;
    vector<int> domainBnd(num_levels);
    for (int i=0; i<num_levels; i++)
        domainBnd[i] = 1;
    dom = MEDDLY::domain::createBottomUp(domainBnd.data(), num_levels);
    for (int i=0; i<num_levels; i++) {
        char buffer[64];
        snprintf(buffer, sizeof(buffer), "x%lu", vorder.lvl2var(i)+1);
        dom->getVar(i + 1)->setName(strdup(buffer));
    }

    // Initialize the MDD forest
    MEDDLY::policies mdd_fp(false);
    mdd_fp.setQuasiReduced();
    mdd_fp.setSparseStorage();
    forestMDD = MEDDLY::forest::create(dom, MEDDLY::SET, MEDDLY::range_type::BOOLEAN, 
                                  MEDDLY::edge_labeling::MULTI_TERMINAL, mdd_fp);

    init_custom_meddly_operators(forestMDD, &pivot_order);
    vzero = zero_vector(forestMDD);
}

/////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge 
sym_canonicalize_gcd(const meddly_context& ctx, MEDDLY::dd_edge A)
{
    // Find the non-trivial divisors of A vectors
    divisors_finder_mdd_op::lut_key lk;
    DIV_FINDER_MDD->computeDDEdge(A, lk);
    const std::vector<int>& divisors = DIV_FINDER_MDD->look_up(lk);

    // divisors appear in descending order
    for (int div : divisors) {
        if (div > 0) {
            MEDDLY::dd_edge U(ctx.forestMDD), D(ctx.forestMDD); // undivisible, divided
            VDIVIDE->computeDDEdge(A, div, U, D);
            A = sym_union(U, D);
        }
    }
    return A;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
sym_s_vectors_at_level(const meddly_context& ctx, const pottier_params_t& pparams,
                       const MEDDLY::dd_edge A, const MEDDLY::dd_edge B, 
                       const size_t level) 
{
    MEDDLY::dd_edge SV(ctx.forestMDD);
    // Get all the domain values at the specified level
    domain_values_enumerator dveA(level), dveB(level);
    dveA.visit(A);
    dveB.visit(B);
    std::vector<int> valA, valB;
    dveA.get_values(valA);
    dveB.get_values(valB);

    // cout << "sym_s_vectors_at_level (level="<<level<<")"<<endl;
    // cout << "valA: "; { for (int i : valA) cout << i << " "; } cout << endl;
    // cout << "valB: "; { for (int j : valB) cout << j << " "; } cout << endl;

    sv_sign svs = pparams.target == compute_target::GRAVER_BASIS ? SVS_UNDECIDED : SVS_POS;
    // svs = SVS_UNDECIDED;
    lv_selector_op* selector_at_level = LV_SELECTOR_OPS->get_op(level);

    // Combine elements from A and B having opposite signs at @level
    for (int i : valA) {
        for (int j : valB) {
            ab_sum_t sum_or_diff;
            int ij_prod = i*j;

            if (ij_prod < 0) {
                sum_or_diff = ab_sum_t::A_PLUS_B; // compute a + b
            }
            else if (ij_prod > 0) {
                // subtract if we are building the graver basis
                if (pparams.target == compute_target::GRAVER_BASIS)
                    sum_or_diff = ab_sum_t::A_MINUS_B; // compute a - b
                else 
                    continue; // do nothing for i and j
            }
            else // ij_prod == 0
                continue;

            // FIXME: ma ij_prod può essere == 0 ????
            if (pparams.target == compute_target::EXTREME_RAYS && 0==ij_prod)
                continue;

            // Select the involved elements
            MEDDLY::dd_edge Ai(ctx.forestMDD), Bj(ctx.forestMDD);
            selector_at_level->computeDDEdge(A, i, Ai);
            selector_at_level->computeDDEdge(B, j, Bj);

            if (pparams.target == compute_target::EXTREME_RAYS) {
                // Find multipliers for i and j
                int mult_Ai = abs(j);
                int mult_Bj = abs(i);
                int g = gcd(mult_Ai, mult_Bj);
                mult_Ai /= g;
                mult_Bj /= g;

                // Hemmecke formula at page 7
                // mult_Bj = mult_Ai / mult_Bj;
                // mult_Ai = 1;
                // if (mult_Bj == 0)
                //     continue;

                // cout << " s-vectors op: "<<i<<"*"<<mult_Ai<<" + "<<j<<"*"<<mult_Bj<<" = "
                //      <<(i*mult_Ai)<<" + "<<(j*mult_Bj)<<" = "
                //      <<(i*mult_Ai + j*mult_Bj) << endl;

                if (mult_Ai)
                    VMULT->computeDDEdge(Ai, mult_Ai, Ai);
                if (mult_Bj)
                    VMULT->computeDDEdge(Bj, mult_Bj, Bj);
            }

            MEDDLY::dd_edge Ai_plus_Bj(ctx.forestMDD);
            S_VECTORS->computeDDEdge(Ai, Bj, false, sum_or_diff, svs, 0, Ai_plus_Bj);
            SV = sym_union(SV, Ai_plus_Bj);
        }
    }

    // if (pparams.target == compute_target::EXTREME_RAYS) {
    //     // Canonicalize the summed entries
    //     SV = sym_canonicalize_gcd(ctx, SV);
    // }
    return SV;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
sym_s_vectors(const meddly_context& ctx, const pottier_params_t& pparams,
              MEDDLY::dd_edge A, const MEDDLY::dd_edge B, 
              const size_t level) 
{
    MEDDLY::dd_edge SV(ctx.forestMDD);
    bool s_vectors_split_at_levels = false;
    if (level != 0 && s_vectors_split_at_levels/*&& pparams.target == compute_target::EXTREME_RAYS*/) {
        // TODO: while s-vectors @ level appears to be faster in benchmark, this should be
        //   disabled through a command-line option.
        // Compute the S-Vectors at the specified level only
        SV = sym_s_vectors_at_level(ctx, pparams, A, B, level);
    }
    else {
        // Compute the S-Vectors over all levels at once

        // a + b
        // When summing vectors, the result will always be in canonical form, as
        // both a and b are. There is therefore no need to decide the sign of the sum,
        // and SVS_POS can be used in all cases.
        S_VECTORS->computeDDEdge(A, B, true, ab_sum_t::A_PLUS_B, SVS_POS, level, SV);

        if (pparams.target == compute_target::GRAVER_BASIS) {
            MEDDLY::dd_edge complSV(ctx.forestMDD);
            // a - b
            // In this case, the sign of the result needs to be checked for canonicity,
            // deciding if the half-basis will store (a-b) or its complement.
            S_VECTORS->computeDDEdge(A, B, true, ab_sum_t::A_MINUS_B, SVS_UNDECIDED, level, complSV);
            SV = sym_union(SV, complSV);
        }
    }
    SV = sym_difference(SV, ctx.vzero);

    // SIGN_CANON_OPS->get_op(true)->computeDDEdge(SV, SV, false);
    if (pparams.target == compute_target::EXTREME_RAYS && pparams.primitive_extremal_rays) {
        // Canonicalize the summed entries (this passes from smallest lattice 
        // representatives to primitive extremal rays)
        SV = sym_canonicalize_gcd(ctx, SV);
    }

    return SV;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
sym_normal_form_extremal_rays(const meddly_context& ctx, const pottier_params_t& pparams,
                              MEDDLY::dd_edge A, const MEDDLY::dd_edge B, 
                              const size_t level) 
{
    // Remove from A the elements having a non-minimal support w.r.t. B
    MEDDLY::dd_edge non_minimal_support(A.getForest());
    // FIXME: using 0 or level doesn't make any difference -> resolved use level
    SUPPORT_INCL_TABLE->get_op(level, true, false)->computeDDEdge(A, B, non_minimal_support, false);
    A = sym_difference(A, non_minimal_support);

    if (!pparams.primitive_extremal_rays) {
        // if extremal rays are not made primitive, we could end up having
        // multiples of the same ray, therefore we need to reduce by subtraction.
        MEDDLY::dd_edge reduced(A.getForest());
        SUPPORT_INCL_TABLE->get_op(level, true, true)->computeDDEdge(non_minimal_support, B, reduced, false);
        A = sym_union(A, reduced);
    }

    // if (!is_emptyset(non_minimal_support)) {
    //     cout << "sym_normal_form_extremal_rays (level="<<level<<")" << endl;
    //     cout << "non_minimal_support:\n" << print_mdd(non_minimal_support, ctx.vorder) << endl;
    //     // cout << "reduced:\n" << print_mdd(reduced, ctx.vorder) << endl;
    // }

    // A = sym_canonicalize_gcd(ctx, A);

    return A;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
sym_normal_form(const meddly_context& ctx, const pottier_params_t& pparams,
                MEDDLY::dd_edge A, const MEDDLY::dd_edge B, 
                const size_t level, bool do_reduction, qnf_op qop) 
{
    // if (dd_cardinality(A) != 1) {
    //     MEDDLY::dd_edge R(ctx.forestMDD);
    //     for (MEDDLY::enumerator path(A); path != 0; ++path) {
    //         const int *ptrs[1];
    //         ptrs[0] = path.getAssignments();
    //         MEDDLY::dd_edge vecA(ctx.forestMDD);
    //         ctx.forestMDD->createEdge(ptrs, 1, vecA);

    //         MEDDLY::dd_edge reduced = sym_normal_form(ctx, pparams, vecA, B, level, do_reduction);
    //         R = sym_union(R, reduced);
    //     }
    //     return R;
    // }

    const size_t normalization_level = (pparams.normalize_by_levels ? level : 0);
    sv_sign svs = pparams.target == compute_target::GRAVER_BASIS ? SVS_UNDECIDED : SVS_POS;
    cmp_sign cs = pparams.target == compute_target::GRAVER_BASIS ? CMP_UNDECIDED : CMP_POS;
    bool is_potentially_eq = (qop == qnf_op::LEQ_NEQ ? true : false);

    if (pparams.target == compute_target::EXTREME_RAYS) {
        return sym_normal_form_extremal_rays(ctx, pparams, A, B, level);
    }

    MEDDLY::dd_edge ALL_I(ctx.forestMDD);

    if (!do_reduction) {
        MEDDLY::dd_edge I(ctx.forestMDD), D(ctx.forestMDD); 
        if (pparams.target == compute_target::GRAVER_BASIS) {
            MEDDLY::dd_edge IP(ctx.forestMDD);
            GET_IRREDUCIBLES->computeDDEdge(A, B, is_potentially_eq, true, svs, CMP_POS, normalization_level, IP, D);
            MEDDLY::dd_edge IN(ctx.forestMDD);
            GET_IRREDUCIBLES->computeDDEdge(A, B, is_potentially_eq, true, svs, CMP_NEG, normalization_level, IN, D);
            I = sym_intersection(IP, IN);
        }
        else {
            GET_IRREDUCIBLES->computeDDEdge(A, B, is_potentially_eq, true, svs, cs, normalization_level, I, D);
        }
        ALL_I = I;
    }
    else {
        // cout << "                     QNF: " << flush;
        while (!is_emptyset(A)) {
            MEDDLY::dd_edge I(ctx.forestMDD), D(ctx.forestMDD);

            if (pparams.target == compute_target::GRAVER_BASIS) {
                MEDDLY::dd_edge IP(ctx.forestMDD), DP(ctx.forestMDD); 
                REDUCE->computeDDEdge(A, B, is_potentially_eq, true, svs, CMP_POS, normalization_level, IP, DP);
                IP = sym_difference(IP, ctx.vzero);
                DP = sym_difference(DP, ctx.vzero);

                MEDDLY::dd_edge IN(ctx.forestMDD), DN(ctx.forestMDD);
                REDUCE->computeDDEdge(A, B, is_potentially_eq, true, svs, CMP_NEG, normalization_level, IN, DN);
                IN = sym_difference(IN, ctx.vzero);
                DN = sym_difference(DN, ctx.vzero);

                I = sym_intersection(IP, IN);
                D = sym_union(DP, DN);
            }
            else {
                REDUCE->computeDDEdge(A, B, is_potentially_eq, true, svs, cs, normalization_level, I, D);
                I = sym_difference(I, ctx.vzero);
                D = sym_difference(D, ctx.vzero);
            }
            // cout << dd_cardinality(D) << " " << flush;



            // cout << "QNF:\n";
            // cout << "A:\n" << print_mdd(A, ctx.vorder) << endl;
            // cout << "B:\n" << print_mdd(B, ctx.vorder) << endl;
            // cout << "I:\n" << print_mdd(I2, ctx.vorder) << endl;
            // cout << "R:\n" << print_mdd(R, ctx.vorder) << endl;
            // cout << "D:\n" << print_mdd(D, ctx.vorder) << endl;

            // if (dd_cardinality(D) > 0) {
            //     cout << "A:\n" << print_mdd(A, ctx.vorder) << endl;
            //     cout << "B:\n" << print_mdd(B, ctx.vorder) << endl;
            //     cout << "I:\n" << print_mdd(I, ctx.vorder) << endl;
            //     cout << "D:\n" << print_mdd(D, ctx.vorder) << endl;
            //     exit(2);
            // }

            ALL_I = sym_union(ALL_I, I);

            // if (!do_reduction)
            //     break;

            // if (pparams.target == compute_target::GRAVER_BASIS)
            //     SIGN_CANON->computeDDEdge(D, D, false); // Not needed as it is done by REDUCE

            A = D;
        }
    }
    // cout << endl;
    return ALL_I;

    // MEDDLY::dd_edge Aprev(A.getForest());
    // do {
    //     Aprev = A;
    //     // Remove from A the elements less_eq_squared B, and add the reduced vectors
    //     MEDDLY::dd_edge reducibles(A.getForest());
    //     // lesseq_sq_op->computeDDEdge(A, B, reducibles, false);
    //     LEQ_NEQ_SQ_COMPARE->computeDDEdge(A, B, true, true, normalization_level, reducibles);
    //     A = sym_difference(A, reducibles);

    //     // if (pparams.half_basis) {
    //     //     // take only the half Graver basis
    //     //     MEDDLY::dd_edge A2(ctx.forestMDD);
    //     //     SIGN_CANON_OPS->get_op(true)->computeDDEdge(A, A2, false);
    //     //     assert(A == A2);
    //     // }

    //     // cout << "reducibles:\n" << print_mdd(reducibles, ctx.vorder) << endl;
    //     if (do_reduction) {
    //         MEDDLY::dd_edge reduced(A.getForest());
    //         LEQ_NEQ_SQ_SUBTRACT->computeDDEdge(A, B, true, true, normalization_level, reduced);
    //         // reduce_sq_op->computeDDEdge(reducibles, B, reduced, false);
    //         // cout << "reduced:\n" << print_mdd(reduced, ctx.vorder) << endl;

    //         A = sym_union(A, reduced);
    //     }
    //     else break; // no need to fix-point
    // } 
    // while (A != Aprev);
    // return A;
}

/////////////////////////////////////////////////////////////////////////////////////////

// Print start of per-iteration banner
static void
pottier_iter_banner_start(const meddly_context& ctx, 
                          const pottier_params_t& pparams,
                          const size_t level,
                          const size_t rem_neg_step,
                          const size_t iter) 
{
    if (pparams.very_verbose) {
        cout << "\n\n";
        sep(pparams);
    }
    if (rem_neg_step != size_t(-1)) {
        bool is_first_level = (level != 0 && level == ctx.pivot_order.var2lvl(0) + 1);
        if (pparams.very_verbose || ((level==0 || is_first_level) && iter==0))
            cout << "Gen=" << left << setw(3) << rem_neg_step;
        else
            cout << "       ";
    }
    if (level != 0) {
        if (pparams.very_verbose || iter == 0) {
            const bool is_pivot = ctx.leading_variables[level-1];
            const std::string& var_name = ctx.forestMDD->getDomain()->getVar(level)->getName();
            cout << "Level=" << level << "["<<var_name<<(is_pivot ? "*" : "")<<"]";
            size_t spaces = 8 - var_name.size();
            if (level >= 1000) spaces--;
            if (level >= 100)  spaces--;
            if (level >= 10)   spaces--;
            if (is_pivot)      spaces--;
            for (size_t i=0; i<spaces; i++)
                cout << " ";
        }
        else
            cout << "                 ";
    }
    cout << "Iter="<< left << setw(3) <<iter;
}

/////////////////////////////////////////////////////////////////////////////////////////

// Pottier algorithm for the computation of the Graver basis in symbolic form
MEDDLY::dd_edge
sym_pottier(const meddly_context& ctx, 
            const pottier_params_t& pparams,
            MEDDLY::dd_edge initGraver, // Graver basis not including N
            MEDDLY::dd_edge N, // new generators
            const size_t level,
            const size_t rem_neg_step)
{
    if (pparams.very_verbose && level>0) {
        const std::string& var_name = ctx.forestMDD->getDomain()->getVar(level)->getName();
        cout << "\n\n\n\n";
        sep(pparams);sep(pparams);
        cout<<"Start of level "<<level<<"["<<var_name<<(ctx.leading_variables[level-1]?"*":"")<<"]"<<endl; 
        cout << "F^init:\n" << print_mdd_lambda(initGraver, ctx.vorder, ctx.pivot_order, level) << endl;
        cout << "N^init:\n" << print_mdd_lambda(N, ctx.vorder, ctx.pivot_order, level) << endl;
        cout << endl;
    }
    if (pparams.graded_order)
        return sym_pottier_grad(ctx, pparams, initGraver, N, level, rem_neg_step);

    bool reduce_C = true;
    degree_type degtype = ((pparams.target == compute_target::EXTREME_RAYS) ? 
                             degree_type::BY_SUPPORT : degree_type::BY_VALUE);

    MEDDLY::dd_edge F(ctx.forestMDD);
    MEDDLY::dd_edge C(ctx.forestMDD);
    MEDDLY::dd_edge prevC(ctx.forestMDD);
    // const MEDDLY::dd_edge empty_set(ctx.forestMDD);
    // sep(pparams);

    // F = N u initG
    F = sym_union(initGraver, N); 

    // cout << "(1) |initGraver|=" << dd_cardinality(initGraver) << ",n="<< initGraver.getNodeCount() << endl;
    // cout << "initGraver:\n" << print_mdd(initGraver, ctx.vorder) << endl;
    // cout << "(1) |N|=" << dd_cardinality(N) << ",n="<< N.getNodeCount() << endl;
    // cout << "N:\n" << print_mdd(N, ctx.vorder) << endl;
    // cout << "(1) |F|=" << dd_cardinality(F) << ",n="<< F.getNodeCount() << endl;
    // cout << "F:\n" << print_mdd(F, ctx.vorder) << endl;
    // assert(dd_cardinality(F) >= std::min(dd_cardinality(initGraver), dd_cardinality(N)));

    if (level != 0 && /*!pparams.normalize_by_levels &&*/ pparams.target!=compute_target::EXTREME_RAYS) {
        // perform the completion procedure to extend to the new column
        // since we are extending the lesseq_sq operation to column j, we need to renormalize F 
        // MEDDLY::dd_edge prevF(ctx.forestMDD);
        F = sym_normal_form(ctx, pparams, F, N, level, true, qnf_op::LEQ_NEQ);

        // MEDDLY::dd_edge removed(ctx.forestMDD), F2(ctx.forestMDD);
        // // COMPL_PROC_OPS->get_op(level)->computeDDEdge(F, F, removed, false);
        // F2 = sym_difference(prevF, F);
        // cout << "initial F:\n" << print_mdd(prevF, ctx.vorder) << endl;
        // if (F != F2) {
        //     cout << "removed:\n" << print_mdd(F2, ctx.vorder) << endl;
        //     cout << "F:\n" << print_mdd(F, ctx.vorder) << endl;
        // }
        // F = F2;
        // cout << "(2) |F|=" << dd_cardinality(F) << ",n="<< F.getNodeCount() << endl;
    }       

    C = sym_s_vectors(ctx, pparams, F, is_emptyset(N) ? F : N, level);
    pparams.perf_C(C);

    size_t iter=0;
    int degree = -1;
    do {
        // Get the subset S of C that will be normalized and furtherly combined
        MEDDLY::dd_edge S(ctx.forestMDD);
        if (pparams.by_degree) { // get the subset of C having the smallest degree
            assert(level != 0);
            int prev_degree = degree;
            SMALLEST_DEGREE_TABLE->get_op(level, degtype)->computeDDEdge(C, degree);
            DEGREE_SELECTOR_TABLE->get_op(level, degtype)->computeDDEdge(C, degree, S);
            // if (!is_emptyset(S)) {
            //     if (!(prev_degree==-1 || prev_degree<=degree)) {
            //         cout << "\n\n\n";
            //         cout << "prev_degree="<<prev_degree<<" degree="<<degree<<endl;
            //         cout << "**C:\n" << print_mdd_lambda(C, ctx.vorder, ctx.pivot_order, level) << endl;
            //         cout << "**S:\n" << print_mdd_lambda(S, ctx.vorder, ctx.pivot_order, level) << endl;
            //     }
            //     assert(prev_degree==-1 || prev_degree<=degree);
            // }
        }
        else { // ignore degrees, take all C at once
            S = C;
        }

        if (pparams.verbose) {
            pottier_iter_banner_start(ctx, pparams, level, rem_neg_step, iter);
            cout << "|F|=" << dd_cardinality(F) << ",n="<< F.getNodeCount();
            cout << "  |C|=" << dd_cardinality(C) << ",n="<< C.getNodeCount();
            if (pparams.by_degree && degree>=0)
                cout << "  |C("<<degree<<")|=" << dd_cardinality(S) << ",n="<< S.getNodeCount();
            cout << endl;
        }

        prevC = C;
        C = sym_difference(C, S);

        // cout << "normal_form" << endl;
        S = sym_normal_form(ctx, pparams, S, F, level, false, qnf_op::LEQ);
        if (!pparams.by_degree)
            S = sym_normal_form(ctx, pparams, S, S, level, false, qnf_op::LEQ_NEQ);
        // S = sym_difference(S, F);

        if (pparams.very_verbose)
            cout << "F:\n" << print_mdd_lambda(F, ctx.vorder, ctx.pivot_order, level) << endl;
        if (is_emptyset(S) && is_emptyset(C))
            break;
        if (pparams.very_verbose) {
            if (!pparams.by_degree) cout << "C:\n";
            else cout << "C("<<degree<<"):\n";
            cout << print_mdd_lambda(S, ctx.vorder, ctx.pivot_order, level) << endl;
        }

        F = sym_union(F, S);
        MEDDLY::dd_edge SV(ctx.forestMDD);
        SV = sym_s_vectors(ctx, pparams, F, S, level);
        pparams.perf_C(SV);

        if (reduce_C) {
            SV = sym_difference(SV, prevC);
            SV = sym_difference(SV, S);
            SV = sym_difference(SV, F);
        }
        // cout << "  |SV|=" << dd_cardinality(SV) << endl;

        // if (pparams.by_degree) { // FIXME: check again this QNF
        //     C = sym_normal_form(ctx, pparams, C, S, level, false);
        // }

        C = sym_union(C, SV);
        // F = sym_union(F, S);

        iter++;
    } while (C != prevC);

    if (pparams.verbose) {
        pottier_iter_banner_start(ctx, pparams, level, rem_neg_step, iter);
        cout << "|F|=" << dd_cardinality(F) << ",n="<< F.getNodeCount();
        if (pparams.by_degree && degree>=0)
            cout << "  max_degree="<<degree<<")";
        cout << endl;
    }

    return F;
}

/////////////////////////////////////////////////////////////////////////////////////////

bool
min_max_degrees(const std::vector<int>& degreesPos,
                const std::vector<int>& degreesNeg,
                bool is_graver,
                int &min_degree, int &max_degree)
{
    min_degree = std::numeric_limits<int>::max();
    max_degree = 0;

    if (is_graver) {
        if (degreesPos.empty() && degreesNeg.empty())
            return false;

        if (!degreesPos.empty()) {
            min_degree = min(min_degree, 2 * degreesPos.back());
            max_degree = max(max_degree, 2 * degreesPos.front());
        }
        if (!degreesNeg.empty()) {
            min_degree = min(min_degree, 2 * degreesNeg.back());
            max_degree = max(max_degree, 2 * degreesNeg.front());
        }
    }
    else {
        if (degreesPos.empty() || degreesNeg.empty())
            return false;

        min_degree = degreesPos.back() + degreesNeg.back();
        max_degree = degreesPos.front() + degreesNeg.front();
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////

// Pottier algorithm for the computation of the Graver basis in symbolic form
// Dynamically generate the S-Vectors by degree, instead of storing them in the C set
MEDDLY::dd_edge
sym_pottier_grad(const meddly_context& ctx, 
                 const pottier_params_t& pparams,
                 MEDDLY::dd_edge initGraver, // Graver basis not including N
                 MEDDLY::dd_edge N, // new generators
                 const size_t level,
                 const size_t rem_neg_step)
{
    assert(level > 0);
    MEDDLY::dd_edge F(ctx.forestMDD);

    // F = N u initG
    F = sym_union(initGraver, N); 

    // TODO: check why EXTREME_RAYS - this should not be there
    if (level != 0 && /*!pparams.normalize_by_levels &&*/ pparams.target!=compute_target::EXTREME_RAYS) {
        // perform the completion procedure to extend to the new column
        // since we are extending the lesseq_sq operation to column j, we need to renormalize F 
        F = sym_normal_form(ctx, pparams, F, N, level, true, qnf_op::LEQ_NEQ);
    }

    if (is_emptyset(F))
        return F;
    degree_type degtype = ((pparams.target == compute_target::EXTREME_RAYS) ? 
                             degree_type::BY_SUPPORT : degree_type::BY_VALUE);

    std::vector<MEDDLY::dd_edge> FposK, FnegK;
    std::vector<bool> havePosK, haveNegK;
    MEDDLY::dd_edge emptySet(ctx.forestMDD);

    int k = -1; // next degree to generate
    int min_gen_degree, max_gen_degree; // range of producible degree
    size_t iter = 0;
    while (true) {
        // rebuild selectors since domain size could change between iterations
        MEDDLY::dd_edge selPos = selector_for_sign_at_level(ctx.forestMDD, level, +1, 1);
        MEDDLY::dd_edge selNeg = selector_for_sign_at_level(ctx.forestMDD, level, -1, 1);
        // separate positives and negatives at @level in F
        MEDDLY::dd_edge Fpos = sym_intersection(F, selPos);
        MEDDLY::dd_edge Fneg = sym_intersection(F, selNeg);
        if (is_emptyset(Fpos) && is_emptyset(Fneg))
            break;

        // Enumerate the degrees of F
        degree_finder_op::lut_key lkP, lkN;
        DEGREE_FINDER_TABLE->get_op(level, degtype)->computeDDEdge(Fpos, lkP);
        DEGREE_FINDER_TABLE->get_op(level, degtype)->computeDDEdge(Fneg, lkN);
        const std::vector<int>& degreesPos = DEGREE_FINDER_TABLE->look_up(lkP);
        const std::vector<int>& degreesNeg = DEGREE_FINDER_TABLE->look_up(lkN);
        // assert(!degreesPos.empty() && !degreesNeg.empty());

        // cout << "|d+|="<<degreesPos.size()<<" ";
        // cout << "|d-|="<<degreesNeg.size();
        // cout << endl;

        if (!min_max_degrees(degreesPos, degreesNeg, pparams.target==compute_target::GRAVER_BASIS,
                             min_gen_degree, max_gen_degree))
            break;

        if (k == -1)
            k = min_gen_degree;
        else if (k > max_gen_degree) 
            break;

        FposK.resize(k+1, emptySet);
        FnegK.resize(k+1, emptySet);
        havePosK.resize(k+1, false);
        haveNegK.resize(k+1, false);

        // if (k == -1) {
        //     // determine smallest producible degree
        //     if (pparams.target == compute_target::GRAVER_BASIS)
        //         k = 2 * min(degreesPos.back(), degreesNeg.back());
        //     else
        //         k = degreesPos.back() + degreesNeg.back();
        // }
        // // always update the maximum producible degree (it changes during the iterations)
        // if (pparams.target == compute_target::GRAVER_BASIS)
        //     max_gen_degree = 2 * max(degreesPos.front(), degreesNeg.front());
        // else
        //     max_gen_degree = degreesPos.front() + degreesNeg.front();

        if (pparams.verbose) {
            pottier_iter_banner_start(ctx, pparams, level, rem_neg_step, iter);
            cout << "|F|=" << dd_cardinality(F) << ",n="<< F.getNodeCount();
            cout << "  k="<<setw(4)<<k;
            cout << "-> "<<setw(4)<<max_gen_degree;
            cout << "|d+|="<<setw(4)<<degreesPos.size();
            if (pparams.very_verbose) {
                for (int d : degreesPos) 
                    cout << d << " ";
            }
            cout << "|d-|="<<setw(4)<<degreesNeg.size();
            if (pparams.very_verbose) {
                for (int d : degreesNeg) 
                    cout << d << " ";
            }
            cout << endl;
        }


        // Generate the S-Vectors of degree k
        // bool added = false;
        bool repeat_for_deg0 = false;
        do {
            bool added_using_deg0 = false;
            for (size_t turn=0; turn<(pparams.target==compute_target::GRAVER_BASIS ? 3 : 1); turn++) {
                const std::vector<int> *degreesI, *degreesJ;
                const MEDDLY::dd_edge *FI, *FJ;
                char signI, signJ;
                ab_sum_t svect_op;
                sv_sign svs;
                std::vector<MEDDLY::dd_edge> *FIK, *FJK;
                std:vector<bool> *haveIK, *haveJK;
                switch (turn) {
                    case 0: // (F+) + (F-)
                        degreesI = &degreesPos;     degreesJ = &degreesNeg;
                        FI = &Fpos;                 FJ = &Fneg;
                        signI = '+';                signJ = '-';
                        svect_op = ab_sum_t::A_PLUS_B;
                        svs = SVS_POS;
                        FIK = &FposK;   haveIK = &havePosK;
                        FJK = &FnegK;   haveJK = &haveNegK;
                        break;

                    case 1: // (F+) - (F+)
                        degreesI = degreesJ = &degreesPos;
                        FI = FJ = &Fpos;
                        signI = signJ = '+';
                        svect_op = ab_sum_t::A_MINUS_B;
                        svs = SVS_UNDECIDED;
                        FIK = &FposK;   haveIK = &havePosK;
                        FJK = &FposK;   haveJK = &havePosK;
                        break;

                    case 2: // (F-) - (F-)
                        degreesI = degreesJ = &degreesNeg;
                        FI = FJ = &Fneg;
                        signI = signJ = '-';
                        svect_op = ab_sum_t::A_MINUS_B;
                        svs = SVS_UNDECIDED;
                        FIK = &FnegK;   haveIK = &haveNegK;
                        FJK = &FnegK;   haveJK = &haveNegK;
                        break;
                }

                // Generate the combinations of i+j = k
                for (int i : *degreesI) {
                    for (int j : *degreesJ) {
                        if (i + j == k) {
                            // skip useless pairs for Graver semi-basis when subtracting from the same sets 
                            if (pparams.target==compute_target::GRAVER_BASIS && svect_op==ab_sum_t::A_MINUS_B && i>j)
                                continue;

                            if (repeat_for_deg0 && (i!=0 && j!=0))
                                continue;

                            // if (i==0 || j==0)
                            //     continue;

                            // MEDDLY::dd_edge Fi(ctx.forestMDD), Fj(ctx.forestMDD);

                            if (!(*haveIK)[i] || i==k) {
                                DEGREE_SELECTOR_TABLE->get_op(level, degtype)->computeDDEdge(*FI, i, (*FIK)[i]);
                                (*haveIK)[i] = true;
                            }
                            if (!(*haveJK)[j] || j==k) {
                                DEGREE_SELECTOR_TABLE->get_op(level, degtype)->computeDDEdge(*FJ, j, (*FJK)[j]);
                                (*haveJK)[j] = true;
                            }
                            MEDDLY::dd_edge Fi = (*FIK)[i];
                            MEDDLY::dd_edge Fj = (*FJK)[j];

                            // DEGREE_SELECTOR_TABLE->get_op(level, degtype)->computeDDEdge(*FI, i, Fi);
                            // DEGREE_SELECTOR_TABLE->get_op(level, degtype)->computeDDEdge(*FJ, j, Fj);
                            // if (is_emptyset(Fi) || is_emptyset(Fj)) // should not happen...
                            //     continue;

                            // S-Vectors of degree k
                            // we can set level=0 since Fi and Fj are already sign-disjoint, so there is no need
                            // to use the sym_s_vectors_at_level() to split values and signs
                            MEDDLY::dd_edge SV(ctx.forestMDD);
                            S_VECTORS->computeDDEdge(Fi, Fj, true, svect_op, svs, level, SV);
                            SV = sym_difference(SV, ctx.vzero);
                            pparams.perf_C(SV);
                            if (pparams.target == compute_target::EXTREME_RAYS && pparams.primitive_extremal_rays) {
                                // Canonicalize the summed entries (this passes from smallest lattice 
                                // representatives to primitive extremal rays)
                                SV = sym_canonicalize_gcd(ctx, SV);
                            }

                            // if (repeat_for_deg0) {
                            //     cout << "rep0 ";
                            //     cout << "  i="<<i<<" |Fi"<<signI<<"|=" << dd_cardinality(Fi) << ",n="<< Fi.getNodeCount();
                            //     cout << "  j="<<j<<" |Fj"<<signJ<<"|=" << dd_cardinality(Fj) << ",n="<< Fj.getNodeCount();
                            //     cout << "  k="<<k<<" |SV|=" << dd_cardinality(SV) << ",n="<< SV.getNodeCount();
                            //     cout << endl;
                            //     cout << "F"<<signI<<"("<<i<<"):\n" << print_mdd_lambda(Fi, ctx.vorder, ctx.pivot_order, level);
                            //     cout << "F"<<signJ<<"("<<j<<"):\n" << print_mdd_lambda(Fj, ctx.vorder, ctx.pivot_order, level);
                            //     cout << "SV("<<k<<"):\n" << print_mdd_lambda(SV, ctx.vorder, ctx.pivot_order, level);
                            //     cout << endl << endl;
                            //     // exit(0);
                            // }

                            if (pparams.very_verbose) {
                                pottier_iter_banner_start(ctx, pparams, level, rem_neg_step, iter);
                                cout << "  i="<<i<<" |Fi"<<signI<<"|=" << dd_cardinality(Fi) << ",n="<< Fi.getNodeCount();
                                cout << "  j="<<j<<" |Fj"<<signJ<<"|=" << dd_cardinality(Fj) << ",n="<< Fj.getNodeCount();
                                cout << "  k="<<k<<" |SV|=" << dd_cardinality(SV) << ",n="<< SV.getNodeCount();
                                cout << endl;
                            }

                            // normalize
                            SV = sym_normal_form(ctx, pparams, SV, F, level, false, qnf_op::LEQ);
                            if (k==0) {
                                SV = sym_normal_form(ctx, pparams, SV, SV, level, false, qnf_op::LEQ_NEQ);
                            }

                            // cout << "  i:"<<i<<" j:"<<j<<"   SV="<<dd_cardinality(SV2)
                            //      << "  normForm(SV)="<<dd_cardinality(SV)<<endl;

                            if (pparams.very_verbose) {
                                cout << "F"<<signI<<"("<<i<<"):\n" << print_mdd_lambda(Fi, ctx.vorder, ctx.pivot_order, level);
                                cout << "F"<<signJ<<"("<<j<<"):\n" << print_mdd_lambda(Fj, ctx.vorder, ctx.pivot_order, level);
                                cout << "SV("<<k<<"):\n" << print_mdd_lambda(SV, ctx.vorder, ctx.pivot_order, level);
                                cout << endl << endl;
                            }

                            if (!is_emptyset(SV)) {
                                MEDDLY::dd_edge oldF(F);
                                F = sym_union(F, SV);
                                if (oldF != F) {
                                    // added = true;
                                    if (i==0 || j==0)
                                        added_using_deg0 = true;
                                }
                            }
                        }
                    }
                }
            }

            if (pparams.very_verbose)
                cout << "F:\n" << print_mdd_lambda(F, ctx.vorder, ctx.pivot_order, level) << endl << endl;

            // repeat if we have not reached the fixed point yet and either k==0, 
            // or for some elemenet k=i+0 (or k=0+j)
            if (added_using_deg0) {
                // repeat for the same k until fixpoint is reached.
                repeat_for_deg0 = true;
                cout << "repeat_for_deg0" << endl;
            }
            else {
                repeat_for_deg0 = false;
            }
            iter++;
        } while (repeat_for_deg0);

        k++;
    }

    if (pparams.verbose && iter==0) {
            pottier_iter_banner_start(ctx, pparams, level, rem_neg_step, iter);
            cout << "|F|=" << dd_cardinality(F) << ",n="<< F.getNodeCount();
            // cout << "  k="<<k<<"("<<min_gen_degree<<"-"<<max_gen_degree<<")";
            
            cout << "  nothing to do." << endl;
        }

    return F;
}

/////////////////////////////////////////////////////////////////////////////////////////

// First outer loop of the Pottier algorithm (project & lift)
MEDDLY::dd_edge
sym_pottier_PnL(const meddly_context& ctx, 
                const pottier_params_t& pparams,
                MEDDLY::dd_edge initGraver, MEDDLY::dd_edge N,
                const std::vector<size_t> *rem_neg_levels, 
                const size_t rem_neg_step)
{
    MEDDLY::dd_edge G(ctx.forestMDD), emptySet(ctx.forestMDD);

    if (pparams.by_levels) {
        
        // Compute by levels, following the pivoting order
        G = initGraver;
        for (size_t var=0; var<ctx.num_levels; var++) {
            // Get the next pivot level that will be completed
            const size_t level = ctx.pivot_order.var2lvl(var) + 1;

            // Get the non-zero initial elements
            MEDDLY::dd_edge nnz_sel = selector_for_nonzeroes_at_level(ctx.forestMDD, level);
            MEDDLY::dd_edge init_level = sym_intersection(N, nnz_sel);
            N = sym_difference(N, init_level);

            if (dd_cardinality(init_level) > 1) {
                cout << "ERROR: selecting more than one initial vector!" << endl;
                const std::string& var_name = ctx.forestMDD->getDomain()->getVar(level)->getName();
                cout << "Before start of level "<<level<<"["<<var_name<<(ctx.leading_variables[level-1]?"*":"")<<"]"<<endl; 
                cout << "init:\n" << print_mdd_lambda(init_level, ctx.vorder, ctx.pivot_order, level) << endl << endl;
                exit(4);
            }

            // if (init_level == emptySet)
            //     init_level = G; // nothing new to add, the initial s-vectors is (G+G)

            // Complete the level
            G = sym_pottier(ctx, pparams, G, init_level, level, rem_neg_step);

            if (pparams.target != compute_target::GRAVER_BASIS) {
                // Remove negative values @level if they are not needed any longer
                if (rem_neg_levels==nullptr || (*rem_neg_levels)[level-1] >= rem_neg_step) {
                    MEDDLY::dd_edge non_neg_sel = selector_for_sign_at_level(ctx.forestMDD, level, +1);
                    G = sym_intersection(G, non_neg_sel);
                }
            }

            // if (var>15) {
            //     exit(6);
            // }
        }

        // // Compute by levels, following the pivoting order
        // G = initGraver;
        // for (size_t var=0; var<ctx.num_levels; var++) {
        //     // Get the next pivot level that will be completed
        //     const size_t level = ctx.pivot_order.var2lvl(var) + 1;

        //     // Complete the level
        //     G = sym_pottier(ctx, pparams, G, N, level, rem_neg_step);
            
        //     if (pparams.target != compute_target::GRAVER_BASIS) {
        //         // Remove negative values @level if they are not needed any longer
        //         if (rem_neg_levels==nullptr || (*rem_neg_levels)[level-1] >= rem_neg_step) {
        //             // // canonicalize vectors that are <0 at @level, and 0 below
        //             // MEDDLY::dd_edge sel = selector_for_sign_at_level(ctx.forestMDD, level, -1, 1);
        //             // for (size_t v=0; v<var; v++) {
        //             //     const size_t lvl = ctx.pivot_order.var2lvl(v) + 1;
        //             //     MEDDLY::dd_edge level_sel = selector_for_value_at_level(ctx.forestMDD, 0, lvl);
        //             //     sel = sym_intersection(sel, level_sel);
        //             // }
        //             // MEDDLY::dd_edge rows_to_canon(ctx.forestMDD);
        //             // rows_to_canon = sym_intersection(G, sel);
        //             // if (!is_emptyset(rows_to_canon)) {
        //             //     cout << " |rows_to_canon|="<<dd_cardinality(rows_to_canon)<<""<<endl;
        //             //     cout << print_mdd(rows_to_canon, ctx.vorder) << endl;
        //             //     G = sym_difference(G, rows_to_canon);
        //             //     VMULT->computeDDEdge(rows_to_canon, -1, rows_to_canon);
        //             //     G = sym_union(G, rows_to_canon);
        //             // }

        //             // cout << "* removing nonzeroes for level "<<level<<endl;
        //             MEDDLY::dd_edge non_neg_sel = selector_for_sign_at_level(ctx.forestMDD, level, +1);
        //             G = sym_intersection(G, non_neg_sel);
        //         }
        //     }

        //     N = G; // N is new only at the first iteration, then combine everything
        //     initGraver = emptySet;
        // }
    }
    else {
        // Compute all levels at once
        G = sym_pottier(ctx, pparams, initGraver, N, 0, rem_neg_step);
    }

    return G;
}

/////////////////////////////////////////////////////////////////////////////////////////

// Second outer loop of the Pottier algorithm (extend and complete)
MEDDLY::dd_edge
sym_pottier_bygen(const meddly_context& ctx, 
                  const pottier_params_t& pparams,
                  const std::vector<std::vector<int>>& lattice_Zgenerators)
{
    MEDDLY::dd_edge G(ctx.forestMDD);
    const MEDDLY::dd_edge empty_set(ctx.forestMDD);
    // We should start with the symmetrized HNF generators only when computing the Graver basis
    const bool make_sign_canonic = (pparams.target == compute_target::GRAVER_BASIS);
    const bool make_gen_sym = false;//(pparams.target != compute_target::GRAVER_BASIS);
    if (lattice_Zgenerators.empty())
        return G;

    if (pparams.by_generators) { // Add one generator at a time, in an outer loop
        // Determine when negative values can be dropped (for Hilbert basis)
        std::vector<size_t> rem_neg_levels = step_for_negative_removal(lattice_Zgenerators);

        std::vector<std::vector<int>> mG(1);
        G = empty_set;
        // Compute by generators, bottom-up w.r.t the integral generators
        for (ssize_t row=lattice_Zgenerators.size()-1; row>=0; --row) {
            mG[0] = lattice_Zgenerators[row];
            // {g} or {g, -g}
            MEDDLY::dd_edge g = mdd_from_vectors(mG, ctx.forestMDD, make_gen_sym);
            if (make_sign_canonic)
                SIGN_CANON->computeDDEdge(g, g, false);

            if (pparams.very_verbose) 
                cout << "\nAdding generator "<<row<<endl;

            // Extend Graver basis G with the new generator g
            if (pparams.graded_EaC)
                G = sym_pottier_EaC_graded(ctx, pparams, G, g, row);
            else
                G = sym_pottier_PnL(ctx, pparams, G, g, &rem_neg_levels, row);

            // For Hilbert basis and extreme rays sets, drop negative vectors (when possible)
            if (pparams.target != compute_target::GRAVER_BASIS /*&& !by_level*/) {
                MEDDLY::dd_edge non_neg_sel = selector_for_nonnegatives(ctx.forestMDD, &rem_neg_levels, row);
                G = sym_intersection(G, non_neg_sel);
            }
            // If extreme rays have not been computed by levels,
            // then we have just the Hilbert basis, from which we can get
            // the extreme rays as the minimal support vectors of the basis.
            if (pparams.target == compute_target::EXTREME_RAYS && !pparams.by_levels) {
                // G = sym_canonicalize_gcd(G, vorder);
                MEDDLY::dd_edge non_minimal_supp(ctx.forestMDD);
                SUPPORT_INCL_TABLE->get_op(0, true, false)->computeDDEdge(G, G, non_minimal_supp, false);
                G = sym_difference(G, non_minimal_supp);
            }

            // Finally, reduce the basis G
            // G = sym_normal_form(ctx, pparams, G, G, 0, false);

    
            // cout << "Gen="<<setw(3)<<row<<" ends with |G|="<<dd_cardinality(G)<<""<<endl;
            // cout << "\nG:\n" << print_mdd(G, ctx.vorder) << "\n\n" << endl;
       }
    }
    else { // Add all generators at once
        // G  or  G u -G
        MEDDLY::dd_edge initF = mdd_from_vectors(lattice_Zgenerators, ctx.forestMDD, make_gen_sym);
        if (make_sign_canonic)
            SIGN_CANON->computeDDEdge(initF, initF, false);

        if (pparams.very_verbose) {
            cout << "Initial Z-generators:\n" << print_mdd_lambda(initF, ctx.vorder, ctx.pivot_order, 0) << endl;
        }

        // Compute all generators at once
        G = sym_pottier_PnL(ctx, pparams, empty_set, initF, nullptr, size_t(-1));
    }

    // For Hilbert basis and extreme rays sets, take only the non-negative vectors
    if (pparams.target != compute_target::GRAVER_BASIS /*&& !by_level*/) {
        MEDDLY::dd_edge non_neg_sel = selector_for_nonnegatives(ctx.forestMDD);//, rem_neg_levels, rem_neg_step);
        G = sym_intersection(G, non_neg_sel);
    }
    // If extreme rays have not been computed by levels,
    // then we have just the Hilbert basis, from which we can get
    // the extreme rays as the minimal support vectors of the basis.
    if (pparams.target == compute_target::EXTREME_RAYS /*&& !by_level*/) {
        // G = sym_canonicalize_gcd(G, vorder);
        MEDDLY::dd_edge non_minimal_supp(ctx.forestMDD);
        SUPPORT_INCL_TABLE->get_op(0, true, false)->computeDDEdge(G, G, non_minimal_supp, false);
        G = sym_difference(G, non_minimal_supp);
    }

    // Finally, apply again the normal form, as the Pottier algorithm only
    // guarantees that G is a superset of the Graver basis, and 
    // we do not proceed by a graded order
    G = sym_normal_form(ctx, pparams, G, G, 0, false, qnf_op::LEQ_NEQ);
    
    return G;
}

//////////////////////////////////////////////////////////////////////////////////////

// Pottier algorithm for the computation of the Graver basis in symbolic form
MEDDLY::dd_edge
sym_pottier_EaC_graded(const meddly_context& ctx, 
                       const pottier_params_t& pparams,
                       MEDDLY::dd_edge initGraver, // Graver basis not including N
                       MEDDLY::dd_edge g, // new generator
                       size_t gen_counter)
{
    if (pparams.verbose) {
        cout<<"Extending with generator g_"<<gen_counter<<endl; 
    }
    if (pparams.very_verbose) {
        sep(pparams);sep(pparams);
        cout << "F^init:\n" << print_mdd(initGraver, ctx.vorder) << endl;
        cout << "g:     \n" << print_mdd(g, ctx.vorder) << endl;
        cout << endl;
    }

    std::vector<MEDDLY::dd_edge> vFd;
    vFd.push_back(initGraver); // F0
    size_t deg = 0, max_deg = 1;

    do {
        if (pparams.very_verbose) {
            cout << "New degree "<<deg << endl;
        }

        deg++;
        MEDDLY::dd_edge Fdeg(ctx.forestMDD);
        if (deg == 1)
            Fdeg = g;
        vFd.push_back(Fdeg);

        // Generate the initial candidates of degree deg
        MEDDLY::dd_edge C(ctx.forestMDD);
        // cout << "   start S-Vectors: "<<flush;
        for (ssize_t kk=1; ssize_t(deg-kk) >= 0; kk++) {
            // cout << kk<< " "<<flush;
            MEDDLY::dd_edge C2(ctx.forestMDD);
            C2 = sym_s_vectors(ctx, pparams, vFd[deg-kk], vFd[kk], 0);
            pparams.perf_C(C2);
            C = sym_union(C, C2);
        }
        // cout<<endl;

        // for (ssize_t kk=1; ssize_t(deg-kk) >= 0; kk++) {
        //     MEDDLY::dd_edge C(ctx.forestMDD);
        //     C = sym_s_vectors(ctx, pparams, vFd[deg-kk], vFd[kk], 0);

        while (!is_emptyset(C)) {
            // cout << "   start normal form: "<<flush;
            for (size_t d=0; d<=min(deg-1, max_deg); d++) {
                // cout << "("<<dd_cardinality(C)<<","<<dd_cardinality(vFd[d])<<") " <<flush;
                C = sym_normal_form(ctx, pparams, C, vFd[d], 0, false, qnf_op::LEQ_NEQ);
                // cout << "* "<<flush;
                C = sym_difference(C, vFd[d]);
            }
            // cout << "D "<<flush;
            C = sym_difference(C, vFd[deg]);
            // cout << "C("<<dd_cardinality(C)<<") " <<flush;
            C = sym_normal_form(ctx, pparams, C, C, 0, false, qnf_op::LEQ_NEQ);
            // cout<<endl;
                
            if (pparams.verbose) {
                cout << "    |F(deg="<<deg<<")|=" << dd_cardinality(vFd[deg]) << ",n="<< vFd[deg].getNodeCount();
                cout << "  |C|=" << dd_cardinality(C) << ",n="<< C.getNodeCount();
                cout << endl;

                if (pparams.very_verbose) {
                    cout << "F(deg="<<deg<<"):\n" << print_mdd(vFd[deg], ctx.vorder) << endl;
                    cout << "C:\n" << print_mdd(C, ctx.vorder) << endl;
                    cout << endl;
                }
            }

            if (!is_emptyset(C)) {
                vFd[deg] = sym_union(vFd[deg], C);
                max_deg = deg;
                C = sym_s_vectors(ctx, pparams, C, initGraver, 0);
                pparams.perf_C(C);
            }
        }
        // }

        if (pparams.verbose) {
                cout << "    |F(deg="<<deg<<")|=" << dd_cardinality(vFd[deg]) << ",n="<< vFd[deg].getNodeCount();
                cout << "  completed."<<endl;
                if (pparams.very_verbose) {
                    cout << "F(deg="<<deg<<"):\n" << print_mdd(vFd[deg], ctx.vorder) << endl;
            }
        }
    }
    while (deg <= 2*max_deg);//(!is_emptyset(vFd[deg]));

    MEDDLY::dd_edge Fout(ctx.forestMDD);
    for (size_t d=0; d<vFd.size(); d++)
        Fout = sym_union(Fout, vFd[d]);

    if (pparams.verbose) {
            cout << "    |G_"<<(gen_counter)<<"|=" << dd_cardinality(Fout) << ",n="<< Fout.getNodeCount();
            cout << "  completed."<<endl;
            if (pparams.very_verbose) {
                cout << "G_"<<(gen_counter)<<":\n" << print_mdd(Fout, ctx.vorder) << endl;
        }
    }

    return Fout;

    // bool reduce_C = true;
    // degree_type degtype = ((pparams.target == compute_target::HILBERT_BASIS) ? 
    //                          degree_type::BY_VALUE : degree_type::BY_SUPPORT);

    // MEDDLY::dd_edge F(ctx.forestMDD);
    // MEDDLY::dd_edge C(ctx.forestMDD);
    // MEDDLY::dd_edge prevC(ctx.forestMDD);
    // // const MEDDLY::dd_edge empty_set(ctx.forestMDD);
    // // sep(pparams);

    // // F = N u initG
    // F = sym_union(initGraver, N); 

    // // cout << "(1) |initGraver|=" << dd_cardinality(initGraver) << ",n="<< initGraver.getNodeCount() << endl;
    // // cout << "initGraver:\n" << print_mdd(initGraver, ctx.vorder) << endl;
    // // cout << "(1) |N|=" << dd_cardinality(N) << ",n="<< N.getNodeCount() << endl;
    // // cout << "N:\n" << print_mdd(N, ctx.vorder) << endl;
    // // cout << "(1) |F|=" << dd_cardinality(F) << ",n="<< F.getNodeCount() << endl;
    // // cout << "F:\n" << print_mdd(F, ctx.vorder) << endl;
    // // assert(dd_cardinality(F) >= std::min(dd_cardinality(initGraver), dd_cardinality(N)));

    // if (level != 0 && !pparams.normalize_by_levels && pparams.target!=compute_target::EXTREME_RAYS) {
    //     // perform the completion procedure to extend to the new column
    //     // since we are extending the lesseq_sq operation to column j, we need to renormalize F 
    //     MEDDLY::dd_edge prevF(ctx.forestMDD);
    //     F = sym_normal_form(ctx, pparams, F, F, level, true);

    //     // MEDDLY::dd_edge removed(ctx.forestMDD), F2(ctx.forestMDD);
    //     // // COMPL_PROC_OPS->get_op(level)->computeDDEdge(F, F, removed, false);
    //     // F2 = sym_difference(prevF, F);
    //     // cout << "initial F:\n" << print_mdd(prevF, ctx.vorder) << endl;
    //     // if (F != F2) {
    //     //     cout << "removed:\n" << print_mdd(F2, ctx.vorder) << endl;
    //     //     cout << "F:\n" << print_mdd(F, ctx.vorder) << endl;
    //     // }
    //     // F = F2;
    // }       

    // C = sym_s_vectors(ctx, pparams, F, is_emptyset(N) ? F : N, level);
    // pparams.perf_C(C);

    // size_t iter=0;
    // int degree = -1;
    // do {
    //     // Get the subset S of C that will be normalized and furtherly combined
    //     MEDDLY::dd_edge S(ctx.forestMDD);
    //     if (pparams.by_degree) { // get the subset of C having the smallest degree
    //         assert(level != 0);
    //         int prev_degree = degree;
    //         SMALLEST_DEGREE_TABLE->get_op(level, degtype)->computeDDEdge(C, degree);
    //         DEGREE_SELECTOR_TABLE->get_op(level, degtype)->computeDDEdge(C, degree, S);
    //         // if (!is_emptyset(S)) {
    //         //     if (!(prev_degree==-1 || prev_degree<=degree)) {
    //         //         cout << "\n\n\n";
    //         //         cout << "prev_degree="<<prev_degree<<" degree="<<degree<<endl;
    //         //         cout << "**C:\n" << print_mdd_lambda(C, ctx.vorder, ctx.pivot_order, level) << endl;
    //         //         cout << "**S:\n" << print_mdd_lambda(S, ctx.vorder, ctx.pivot_order, level) << endl;
    //         //     }
    //         //     assert(prev_degree==-1 || prev_degree<=degree);
    //         // }
    //     }
    //     else { // ignore degrees, take all C at once
    //         S = C;
    //     }

    //     if (pparams.verbose) {
    //         pottier_iter_banner_start(ctx, pparams, level, rem_neg_step, iter);
    //         cout << "|F|=" << dd_cardinality(F) << ",n="<< F.getNodeCount();
    //         cout << "  |C|=" << dd_cardinality(C) << ",n="<< C.getNodeCount();
    //         if (pparams.by_degree && degree>=0)
    //             cout << "  |C("<<degree<<")|=" << dd_cardinality(S) << ",n="<< S.getNodeCount();
    //         cout << endl;
    //     }

    //     prevC = C;
    //     C = sym_difference(C, S);

    //     // cout << "normal_form" << endl;
    //     if (!pparams.by_degree)
    //         S = sym_normal_form(ctx, pparams, S, S, level, true);
    //     S = sym_normal_form(ctx, pparams, S, F, level, true);
    //     S = sym_difference(S, F);

    //     if (pparams.very_verbose)
    //         cout << "F:\n" << print_mdd_lambda(F, ctx.vorder, ctx.pivot_order, level) << endl;
    //     if (is_emptyset(S) && is_emptyset(C))
    //         break;
    //     if (pparams.very_verbose) {
    //         if (!pparams.by_degree) cout << "C:\n";
    //         else cout << "C("<<degree<<"):\n";
    //         cout << print_mdd_lambda(S, ctx.vorder, ctx.pivot_order, level) << endl;
    //     }

    //     F = sym_union(F, S);
    //     MEDDLY::dd_edge SV(ctx.forestMDD);
    //     SV = sym_s_vectors(ctx, pparams, F, S, level);
    //     pparams.perf_C(SV);

    //     if (reduce_C) {
    //         SV = sym_difference(SV, prevC);
    //         SV = sym_difference(SV, S);
    //         SV = sym_difference(SV, F);
    //     }
    //     // cout << "  |SV|=" << dd_cardinality(SV) << endl;

    //     // if (pparams.by_degree) { // FIXME: check again this QNF
    //     //     C = sym_normal_form(ctx, pparams, C, S, level, false);
    //     // }

    //     C = sym_union(C, SV);
    //     // F = sym_union(F, S);

    //     iter++;
    // } while (C != prevC);

    // return F;
}

//////////////////////////////////////////////////////////////////////////////////////














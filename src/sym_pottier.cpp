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
    dom = MEDDLY::createDomainBottomUp(domainBnd.data(), num_levels);
    for (int i=0; i<num_levels; i++) {
        char buffer[64];
        snprintf(buffer, sizeof(buffer), "x%lu", vorder.lvl2var(i)+1);
        dom->useVar(i + 1)->setName(strdup(buffer));
    }

    // Initialize the MDD forest
    MEDDLY::policies mdd_fp(false);
    mdd_fp.setQuasiReduced();
    mdd_fp.setSparseStorage();
    forestMDD = dom->createForest(false, MEDDLY::range_type::BOOLEAN, 
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
    if (level != 0 /*&& pparams.target == compute_target::EXTREME_RAYS*/) {
        // TODO: while s-vectors @ level appears to be faster in benchmark, this should be
        //   disabled through a command-line option.
        // Compute the S-Vectors at the specified level only
        SV = sym_s_vectors_at_level(ctx, pparams, A, B, level);
    }
    else {
        // Compute the S-Vectors over all levels at once
        sv_sign svs = pparams.target == compute_target::GRAVER_BASIS ? SVS_UNDECIDED : SVS_POS;
        // a + b
        // TODO: no need for svs even in Graver, use SVS_POS always.
        S_VECTORS->computeDDEdge(A, B, true, ab_sum_t::A_PLUS_B, svs, level, SV);

        if (pparams.target == compute_target::GRAVER_BASIS) {
            MEDDLY::dd_edge complSV(ctx.forestMDD);
            // a - b
            S_VECTORS->computeDDEdge(A, B, true, ab_sum_t::A_MINUS_B, svs, level, complSV);
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
                const size_t level, bool do_reduction) 
{
    const size_t normalization_level = (pparams.normalize_by_levels ? level : 0);
    sv_sign svs = pparams.target == compute_target::GRAVER_BASIS ? SVS_UNDECIDED : SVS_POS;
    cmp_sign cs = pparams.target == compute_target::GRAVER_BASIS ? CMP_UNDECIDED : CMP_POS;

    if (pparams.target == compute_target::EXTREME_RAYS) {
        return sym_normal_form_extremal_rays(ctx, pparams, A, B, level);
    }

    MEDDLY::dd_edge I(ctx.forestMDD);
    while (!is_emptyset(A)) {
        MEDDLY::dd_edge I2(ctx.forestMDD), R(ctx.forestMDD), D(ctx.forestMDD);
        REDUCE->computeDDEdge(A, B, true, true, svs, cs, normalization_level, R, D);
        I2 = sym_difference(A, R);

        // FIXME: REDUCE può tornare solo R e D, non serve I

        // cout << "QNF:\n";
        // cout << "A:\n" << print_mdd(A, ctx.vorder) << endl;
        // cout << "B:\n" << print_mdd(B, ctx.vorder) << endl;
        // cout << "I:\n" << print_mdd(I2, ctx.vorder) << endl;
        // cout << "R:\n" << print_mdd(R, ctx.vorder) << endl;
        // cout << "D:\n" << print_mdd(D, ctx.vorder) << endl;

        I = sym_union(I, I2);
        // TODO: check I does not contain the zero vector.

        // if (pparams.target == compute_target::GRAVER_BASIS)
        //     SIGN_CANON->computeDDEdge(D, D, false); // Not needed as it is done by REDUCE

        A = D;
    }
    return I;

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
            const char* var_name = ctx.forestMDD->getDomain()->getVar(level)->getName();
            cout << "Level=" << level << "["<<var_name<<"]";
            size_t spaces = 10 - strlen(var_name);
            if (level < 1000) spaces--;
            if (level < 100)  spaces--;
            if (level < 10)   spaces--;
            for (size_t i=0; i<spaces; i++)
                cout << " ";
        }
        else
            cout << "                ";
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
    if (pparams.dynamic_svectors)
        return sym_pottier_grad(ctx, pparams, initGraver, N, level, rem_neg_step);

    bool reduce_C = true;
    degree_type degtype = ((pparams.target == compute_target::HILBERT_BASIS) ? 
                             degree_type::BY_VALUE : degree_type::BY_SUPPORT);

    MEDDLY::dd_edge F(ctx.forestMDD);
    MEDDLY::dd_edge C(ctx.forestMDD);
    MEDDLY::dd_edge prevC(ctx.forestMDD);
    // const MEDDLY::dd_edge empty_set(ctx.forestMDD);
    // sep(pparams);

    // F = N u initG
    F = sym_union(initGraver, N); 
    // if (pparams.by_generators) {
    //     F = initGraver;
    // }
    // else {
    //     // F = N u initG
    //     F = sym_union(initGraver, N); 
    // }

    // cout << "(1) |initGraver|=" << initGraver.getCardinality() << ",n="<< initGraver.getNodeCount() << endl;
    // cout << "initGraver:\n" << print_mdd(initGraver, ctx.vorder) << endl;
    // cout << "(1) |N|=" << N.getCardinality() << ",n="<< N.getNodeCount() << endl;
    // cout << "N:\n" << print_mdd(N, ctx.vorder) << endl;
    // cout << "(1) |F|=" << F.getCardinality() << ",n="<< F.getNodeCount() << endl;
    // cout << "F:\n" << print_mdd(F, ctx.vorder) << endl;
    // assert(F.getCardinality() >= std::min(initGraver.getCardinality(), N.getCardinality()));

    if (level != 0 && !pparams.normalize_by_levels && pparams.target!=compute_target::EXTREME_RAYS) {
        // perform the completion procedure to extend to the new column
        // since we are extending the lesseq_sq operation to column j, we need to renormalize F 
        MEDDLY::dd_edge prevF(ctx.forestMDD);
        // do {
        //     prevF = F;
            F = sym_normal_form(ctx, pparams, F, F, level, false);
        // } while (F != prevF);

        // MEDDLY::dd_edge removed(ctx.forestMDD), F2(ctx.forestMDD);
        // // COMPL_PROC_OPS->get_op(level)->computeDDEdge(F, F, removed, false);
        // F2 = sym_difference(prevF, F);
        // cout << "initial F:\n" << print_mdd(prevF, ctx.vorder) << endl;
        // if (F != F2) {
        //     cout << "removed:\n" << print_mdd(F2, ctx.vorder) << endl;
        //     cout << "F:\n" << print_mdd(F, ctx.vorder) << endl;
        // }
        // F = F2;
    }       

    // TODO: se si va per generatori, avremmo potuto non fare subito la F=union(initG, N), per non ricombinare con N
    C = sym_s_vectors(ctx, pparams, F, N, level);
    pparams.perf_C(C);

    // if (pparams.by_generators) {
    //     F = sym_union(F, N);
    // }

    size_t iter=0;
    do {
        // Get the subset S of C that will be normalized and furtherly combined
        MEDDLY::dd_edge S(ctx.forestMDD);
        int degree = -1;
        if (pparams.by_degree) { // get the subset of C having the smallest degree
            assert(level != 0);
            SMALLEST_DEGREE_TABLE->get_op(level, degtype)->computeDDEdge(C, degree);
            DEGREE_SELECTOR_TABLE->get_op(level, degtype)->computeDDEdge(C, degree, S);
        }
        else { // ignore degrees, take all C at once
            S = C;
        }

        if (pparams.verbose) {
            pottier_iter_banner_start(ctx, pparams, level, rem_neg_step, iter);
            cout << "|F|=" << F.getCardinality() << ",n="<< F.getNodeCount();
            cout << "  |C|=" << C.getCardinality() << ",n="<< C.getNodeCount();
            if (pparams.by_degree && degree>=0)
                cout << "  |C("<<degree<<")|=" << S.getCardinality() << ",n="<< S.getNodeCount();
            cout << endl;
        }

        prevC = C;
        C = sym_difference(C, S);


        // cout << "normal_form" << endl;
        if (!pparams.by_degree)
            S = sym_normal_form(ctx, pparams, S, S, level, false);
        S = sym_normal_form(ctx, pparams, S, F, level, false);
        S = sym_difference(S, F);

        if (pparams.very_verbose)
            cout << "F:\n" << print_mdd(F, ctx.vorder) << endl;
        if (is_emptyset(S) && is_emptyset(C))
            break;
        if (pparams.very_verbose) {
            cout << "C:\n" << print_mdd(S, ctx.vorder) << endl;
        }

        MEDDLY::dd_edge SV(ctx.forestMDD);
        SV = sym_s_vectors(ctx, pparams, F, S, level);
        pparams.perf_C(SV);
        if (reduce_C) {
            SV = sym_difference(SV, prevC);
            SV = sym_difference(SV, S);
            SV = sym_difference(SV, F);
        }
        // cout << "  |SV|=" << SV.getCardinality() << endl;

        // if (pparams.by_degree) { // FIXME: check again this QNF
        //     C = sym_normal_form(ctx, pparams, C, S, level, false);
        // }

        C = sym_union(C, SV);
        F = sym_union(F, S);

        iter++;
    } while (C != prevC);

    return F;
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

    if (pparams.very_verbose) {
        cout << "\n\n\n\nStart of level "<<level<<endl; 
        cout << "F:\n" << print_mdd(F, ctx.vorder) << endl;
        cout << endl;
    }

    if (level != 0 && !pparams.normalize_by_levels && pparams.target!=compute_target::EXTREME_RAYS) {
        // perform the completion procedure to extend to the new column
        // since we are extending the lesseq_sq operation to column j, we need to renormalize F 
        F = sym_normal_form(ctx, pparams, F, F, level, false);
    }

    if (is_emptyset(F))
        return F;
    degree_type degtype = ((pparams.target == compute_target::HILBERT_BASIS) ? 
                             degree_type::BY_VALUE : degree_type::BY_SUPPORT);

    int k = -1; // next degree to generate
    size_t iter = 0;
    while (true) {
        // rebuild selectors since domain size could change between iterations
        MEDDLY::dd_edge selPos = selector_for_sign_at_level(ctx.forestMDD, level, +1, 1);
        MEDDLY::dd_edge selNeg = selector_for_sign_at_level(ctx.forestMDD, level, -1, 1);
        // separate positives and negatives at @level in F
        MEDDLY::dd_edge Fpos = sym_intersection(F, selPos);
        MEDDLY::dd_edge Fneg = sym_intersection(F, selNeg);
        if (is_emptyset(Fpos) || is_emptyset(Fneg))
            break;

        // Enumerate the degrees of F
        degree_finder_op::lut_key lkP, lkN;
        DEGREE_FINDER_TABLE->get_op(level, degtype)->computeDDEdge(Fpos, lkP);
        DEGREE_FINDER_TABLE->get_op(level, degtype)->computeDDEdge(Fneg, lkN);
        const std::vector<int>& degreesPos = DEGREE_FINDER_TABLE->look_up(lkP);
        const std::vector<int>& degreesNeg = DEGREE_FINDER_TABLE->look_up(lkN);
        assert(!degreesPos.empty() && !degreesNeg.empty());

        if (k == -1)
            k = degreesPos.back() + degreesNeg.back(); // smallest producible degree
        if (k > degreesPos.front() + degreesNeg.front()) // largest producible degree
            break;

        if (pparams.verbose) {
            pottier_iter_banner_start(ctx, pparams, level, rem_neg_step, iter);
            cout << "|F|=" << F.getCardinality() << ",n="<< F.getNodeCount();
            cout << "  k="<<setw(4)<<k;
            cout << "-> "<<setw(4)<<(degreesPos.front() + degreesNeg.front());
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
        bool added = false;
        for (int i : degreesPos) {
            for (int j : degreesNeg) {
                if (i + j == k) {
                    MEDDLY::dd_edge Fi(ctx.forestMDD), Fj(ctx.forestMDD), SV(ctx.forestMDD);

                    DEGREE_SELECTOR_TABLE->get_op(level, degtype)->computeDDEdge(Fpos, i, Fi);
                    if (is_emptyset(Fi)) 
                        continue;
                    DEGREE_SELECTOR_TABLE->get_op(level, degtype)->computeDDEdge(Fneg, j, Fj);
                    if (is_emptyset(Fj)) 
                        continue;

                    // S-Vectors of degree k
                    // we can set level=0 since Fi and Fj are already sign-disjoint, so there is no need
                    // to use the sym_s_vectors_at_level() to split values and signs
                    SV = sym_s_vectors(ctx, pparams, Fi, Fj, 
                                /*level*/ pparams.target == compute_target::HILBERT_BASIS ? 0 : level);
                    pparams.perf_C(SV);
                    // normalize
                    SV = sym_normal_form(ctx, pparams, SV, F, level, false);
                    // cout << "  i:"<<i<<" j:"<<j<<"   SV="<<SV2.getCardinality()
                    //      << "  normForm(SV)="<<SV.getCardinality()<<endl;

                    if (pparams.very_verbose) {
                        cout << "F+("<<i<<"):\n" << print_mdd(Fi, ctx.vorder);
                        cout << "F-("<<j<<"):\n" << print_mdd(Fj, ctx.vorder);
                        cout << "SV("<<k<<"):\n" << print_mdd(SV, ctx.vorder);
                        cout << endl << endl;
                    }

                    if (!is_emptyset(SV)) {
                        MEDDLY::dd_edge oldF(F);
                        F = sym_union(F, SV);
                        if (oldF != F)
                            added = true;
                    }
                }
            }
        }

        if (pparams.very_verbose)
            cout << "F:\n" << print_mdd(F, ctx.vorder) << endl << endl;

        if (added && k==0) {
            // repeat for k=0
            // FIXME: this may happen also when k>0, but there are
            // 0-degree vectors for i and j.
        }
        else
            k++;
        iter++;
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

            if (init_level == emptySet)
                init_level = G; // nothing new to add, the initial s-vectors is (G+G)

            // Complete the level
            G = sym_pottier(ctx, pparams, G, init_level, level, rem_neg_step);

            if (pparams.target != compute_target::GRAVER_BASIS) {
                // Remove negative values @level if they are not needed any longer
                if (rem_neg_levels==nullptr || (*rem_neg_levels)[level-1] >= rem_neg_step) {
                    MEDDLY::dd_edge non_neg_sel = selector_for_sign_at_level(ctx.forestMDD, level, +1);
                    G = sym_intersection(G, non_neg_sel);
                }
            }
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
        //             //     cout << " |rows_to_canon|="<<rows_to_canon.getCardinality()<<""<<endl;
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
    const bool make_gen_sym = (pparams.target != compute_target::GRAVER_BASIS);
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
    
            // cout << "Gen="<<setw(3)<<row<<" ends with |G|="<<G.getCardinality()<<""<<endl;
            // cout << "\nG:\n" << print_mdd(G, ctx.vorder) << "\n\n" << endl;
       }
    }
    else { // Add all generators at once
        // G  or  G u -G
        MEDDLY::dd_edge initF = mdd_from_vectors(lattice_Zgenerators, ctx.forestMDD, make_gen_sym);
        if (make_sign_canonic)
            SIGN_CANON->computeDDEdge(initF, initF, false);
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
    G = sym_normal_form(ctx, pparams, G, G, 0, false);
    
    return G;
}

//////////////////////////////////////////////////////////////////////////////////////














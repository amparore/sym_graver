/////////////////////////////////////////////////////////////////////////////////////////
#ifndef __SYMBOLIC_POTTIER_H__
#define __SYMBOLIC_POTTIER_H__
/////////////////////////////////////////////////////////////////////////////////////////

#include "dd_operations.h"
#include "matrix.h"
#include "variable_order.h"

/////////////////////////////////////////////////////////////////////////////////////////

struct perf_counters_t {
    size_t counter_C = 0;
};

/////////////////////////////////////////////////////////////////////////////////////////

enum class compute_target {
    GRAVER_BASIS, HILBERT_BASIS, EXTREME_RAYS
};

struct pottier_params_t {
    bool by_levels = true;
    bool normalize_by_levels = false;
    bool by_generators = false;
    bool verbose = false;
    bool very_verbose = false;
    bool verbose_for_basis = false;
    bool by_degree = true;
    bool dynamic_svectors = false;
    compute_target target = compute_target::HILBERT_BASIS;

    // At which step in the by-gen process it becomes possible 
    // to drop all negatives at a given level (Hilbert/exray only)
    std::vector<size_t> rem_neg_levels;

    perf_counters_t *perf = nullptr;

    inline void perf_C(const MEDDLY::dd_edge& C) const {
        if (perf) {  perf->counter_C += C.getCardinality();  }
    }

    inline bool verbose_show_mat(const std::vector<std::vector<int>>& A) const {
        size_t n = A.size();
        if (n > 0) 
            n = std::max(n, A[0].size());
        return (verbose && n<50) || very_verbose;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////

// keep together the basic state variables needed to use meddly
struct meddly_context {
    const variable_order  vorder;       // mapping of the original variables -> MDD levels
    const variable_order  pivot_order;  // order of the pivot completion when  going by level
    MEDDLY::forest       *forestMDD = nullptr; // MDD forest
    const size_t          num_levels;   // number of variables/levels
    MEDDLY::dd_edge       vzero;        // zero vector


    meddly_context(size_t num_levels, const variable_order& vorder, const variable_order& pivot_order)
    /**/ : num_levels(num_levels), vorder(vorder), pivot_order(pivot_order) {}

    void initialize(size_t meddly_cache_size);
};

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge 
sym_union(const MEDDLY::dd_edge& A, const MEDDLY::dd_edge& B);

MEDDLY::dd_edge 
sym_intersection(const MEDDLY::dd_edge& A, const MEDDLY::dd_edge& B);

MEDDLY::dd_edge 
sym_difference(const MEDDLY::dd_edge& A, const MEDDLY::dd_edge& B);

/////////////////////////////////////////////////////////////////////////////////////////

// print separator
void sep(const pottier_params_t& pparams);

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge 
sym_canonicalize_gcd(const meddly_context& ctx, MEDDLY::dd_edge A);

MEDDLY::dd_edge
sym_s_vectors_at_level(const meddly_context& ctx, const pottier_params_t& pparams,
                       const MEDDLY::dd_edge A, const MEDDLY::dd_edge B, 
                       const size_t level);

MEDDLY::dd_edge
sym_s_vectors(const meddly_context& ctx, const pottier_params_t& pparams,
              MEDDLY::dd_edge A, const MEDDLY::dd_edge B, 
              const size_t level);

MEDDLY::dd_edge
sym_normal_form_extremal_rays(const meddly_context& ctx, const pottier_params_t& pparams,
                              MEDDLY::dd_edge A, const MEDDLY::dd_edge B, 
                              const size_t level);

MEDDLY::dd_edge
sym_normal_form(const meddly_context& ctx, const pottier_params_t& pparams,
                MEDDLY::dd_edge A, const MEDDLY::dd_edge B, 
                const size_t level, bool do_reduction);

// Pottier algorithm for the computation of the Graver basis in symbolic form
MEDDLY::dd_edge
sym_pottier(const meddly_context& ctx, 
            const pottier_params_t& pparams,
            MEDDLY::dd_edge initGraver, // Graver basis not including N
            MEDDLY::dd_edge N, // new generators
            const size_t level, // step when going by levels, or 0
            const size_t rem_neg_step); // step when going by gens, or -1

// Pottier algorithm for the computation of the Graver basis in symbolic form
// Dynamically generate the S-Vectors by degree, istead of storing them in the C set
MEDDLY::dd_edge
sym_pottier_grad(const meddly_context& ctx, 
                 const pottier_params_t& pparams,
                 MEDDLY::dd_edge initGraver, // Graver basis not including N
                 MEDDLY::dd_edge N, // new generators
                 const size_t level,
                 const size_t rem_neg_step);

// First outer loop of the Pottier algorithm (project & lift)
MEDDLY::dd_edge
sym_pottier_PnL(const meddly_context& ctx, 
                const pottier_params_t& pparams,
                MEDDLY::dd_edge initGraver, MEDDLY::dd_edge N,
                const std::vector<size_t> *rem_neg_levels, 
                const size_t rem_neg_step);

// Second outer loop of the Pottier algorithm (extend and complete)
MEDDLY::dd_edge
sym_pottier_bygen(const meddly_context& ctx, 
                  const pottier_params_t& pparams,
                  const std::vector<std::vector<int>>& lattice_Zgenerators);

/////////////////////////////////////////////////////////////////////////////////////////
#endif // __SYMBOLIC_POTTIER_H__
/////////////////////////////////////////////////////////////////////////////////////////

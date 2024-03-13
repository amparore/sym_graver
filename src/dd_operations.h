/////////////////////////////////////////////////////////////////////////////////////////
#ifndef __OPERATORS_H__
#define __OPERATORS_H__
/////////////////////////////////////////////////////////////////////////////////////////
#include <meddly/meddly.h>
#include <meddly/ct_entry_key.h>
#include <meddly/ct_entry_result.h>
#include <meddly/compute_table.h>
/////////////////////////////////////////////////////////////////////////////////////////

#ifdef DD_OPERATIONS_IMPL
# define DD_EXTERN
#else
# define DD_EXTERN extern
#endif

/////////////////////////////////////////////////////////////////////////////////////////
// Custom MEDDLY operators
/////////////////////////////////////////////////////////////////////////////////////////

class variable_order;
// class s_vectors_opname;
// class s_vectors_table;
// class s_vectors;

// DD_EXTERN s_vectors_opname *S_VECTORS_OPNAME;
// DD_EXTERN s_vectors_table *S_VECTORS_OPS;

class s_vectors_opname;
class s_vectors;

DD_EXTERN s_vectors_opname *S_VECTORS_OPNAME;
DD_EXTERN s_vectors *S_VECTORS;


// class lesseq_sq_opname;
// class lesseq_sq_table;
// class lesseq_sq;

// DD_EXTERN lesseq_sq_opname *LESSEQ_SQ_OPNAME;
// DD_EXTERN lesseq_sq_table *LESSEQ_SQ_OPS;

class leq_neq_sq;
class leq_neq_sq_opname;

DD_EXTERN leq_neq_sq_opname *LEQ_NEQ_SQ_OPNAME;
DD_EXTERN leq_neq_sq *LEQ_NEQ_SQ_COMPARE;
DD_EXTERN leq_neq_sq *LEQ_NEQ_SQ_SUBTRACT;


class reduce_opname;
class reduce;
DD_EXTERN reduce_opname *REDUCE_OPNAME;
DD_EXTERN reduce *REDUCE;
DD_EXTERN reduce *GET_IRREDUCIBLES;

// class reduce3_opname;
// class reduce3;
// DD_EXTERN reduce3_opname *REDUCE3_OPNAME;
// DD_EXTERN reduce3 *REDUCE3;

// class compl_proc_opname;
// class compl_proc_table;
// class compl_proc;

// DD_EXTERN compl_proc_opname *COMPL_PROC_OPNAME;
// DD_EXTERN compl_proc_table *COMPL_PROC_OPS;

class sign_canon_mdd_opname;
class sign_canon_mdd_op;
DD_EXTERN sign_canon_mdd_opname* SIGN_CANON_OPNAME;
DD_EXTERN sign_canon_mdd_op* SIGN_CANON;

class vmult_opname;
class vmult_op;
DD_EXTERN vmult_opname* VMULT_OPNAME;
DD_EXTERN vmult_op* VMULT;

class vdivide_mdd_op;
class vdivide_mdd_opname;
DD_EXTERN vdivide_mdd_opname* VDIVIDE_OPNAME;
DD_EXTERN vdivide_mdd_op* VDIVIDE;

class divisors_finder_mdd_op;
class divisors_finder_mdd_opname;

DD_EXTERN divisors_finder_mdd_opname* DIV_FINDER_MDD_OPNAME;
DD_EXTERN divisors_finder_mdd_op* DIV_FINDER_MDD;

class support_inclusion_opname;
class support_inclusion_table;
DD_EXTERN support_inclusion_opname *SUPPORT_INCL_OPNAME;
DD_EXTERN support_inclusion_table *SUPPORT_INCL_TABLE;

class lv_selector_opname;
class lv_selector_table;
DD_EXTERN lv_selector_opname* LV_SELECTOR_OPNAME;
DD_EXTERN lv_selector_table* LV_SELECTOR_OPS;

class smallest_degree_opname;
class smallest_degree_table;
DD_EXTERN smallest_degree_opname *SMALLEST_DEGREE_OPNAME;
DD_EXTERN smallest_degree_table *SMALLEST_DEGREE_TABLE;

class degree_finder_opname;
class degree_finder_table;
DD_EXTERN degree_finder_opname *DEGREE_FINDER_OPNAME;
DD_EXTERN degree_finder_table *DEGREE_FINDER_TABLE;

class degree_selector_opname;
class degree_selector_table;
DD_EXTERN degree_selector_opname *DEGREE_SELECTOR_OPNAME;
DD_EXTERN degree_selector_table *DEGREE_SELECTOR_TABLE;

void init_custom_meddly_operators(MEDDLY::forest* forestMDD, const variable_order *pivot_order);

/////////////////////////////////////////////////////////////////////////////////////////
// Encoding of integers using Meddly nodes (restricted to naturals)
/////////////////////////////////////////////////////////////////////////////////////////

inline unsigned ZtoNode(int z) {
    return (z >= 0) ? unsigned(z<<1) : unsigned((-z)<<1)-1;
}
inline int NodeToZ(unsigned n) {
    return (n&1)==0 ? int(n>>1) : -(int((n+1)>>1));
}

/////////////////////////////////////////////////////////////////////////////////////////
class variable_order;

struct mdd_printer {
    MEDDLY::dd_edge dd;
    const variable_order& vorder;
    bool write_header;
    bool write_counts;
    // only for projection problems
    const variable_order* pivot_order;
    const size_t lambda;
};
std::ostream& operator<< (std::ostream&, const mdd_printer&);

inline mdd_printer print_mdd(MEDDLY::dd_edge dd, const variable_order& vorder,
                             bool write_header=true, bool write_counts=false) {
    return mdd_printer{ .dd=dd, .vorder=vorder, .write_header=write_header, 
                        .write_counts=write_counts, .pivot_order=nullptr, .lambda=0 };
}

inline mdd_printer print_mdd_lambda(MEDDLY::dd_edge dd, const variable_order& vorder,
                                    const variable_order& pivot_order, size_t lambda,
                                    bool write_header=true, bool write_counts=false) {
    return mdd_printer{ .dd=dd, .vorder=vorder, .write_header=write_header, 
                        .write_counts=write_counts, .pivot_order=&pivot_order, .lambda=lambda };
}

/////////////////////////////////////////////////////////////////////////////////////////

void check_level_bound(MEDDLY::forest* forest, int level, int bound);

MEDDLY::dd_edge 
mdd_from_vectors(const std::vector<std::vector<int>>& vecs, 
                 MEDDLY::forest* forestMDD, bool make_symmetric);

inline bool is_emptyset(const MEDDLY::dd_edge& e) {
    return e.getNode() == e.getForest()->handleForValue(false);
}

MEDDLY::dd_edge
selector_for_value_at_level(MEDDLY::forest *forestMDD, int value, size_t level);

MEDDLY::dd_edge
selector_for_nonzeroes_at_level(MEDDLY::forest *forestMDD, size_t level);

MEDDLY::dd_edge
selector_for_sign_at_level(MEDDLY::forest *forestMDD, size_t level, int sign, int threshold=0);

MEDDLY::dd_edge
selector_for_nonnegatives(MEDDLY::forest *forestMDD);

MEDDLY::dd_edge
selector_for_nonnegatives(MEDDLY::forest *forestMDD, 
                          const std::vector<size_t> *rmnnz_levels, size_t rmnnz_step);

MEDDLY::dd_edge
selector_for_zeros_up_to_level(MEDDLY::forest *forestMDD, size_t max_level);

MEDDLY::dd_edge
zero_vector(MEDDLY::forest *forestMDD);

void unionNodes(MEDDLY::unpacked_node* C, MEDDLY::node_handle node, 
                unsigned pos, MEDDLY::forest* resF,  
                MEDDLY::binary_operation* mddUnion);

double dd_cardinality(const MEDDLY::dd_edge& dd);

/////////////////////////////////////////////////////////////////////////////////////////
// Save DD encoding Integers to dot/pdf format
void write_dd_as_pdf(const MEDDLY::dd_edge& e, 
                     const char* base_name, bool level_labels, 
                     bool write_terminals);

/////////////////////////////////////////////////////////////////////////////////////////









/////////////////////////////////////////////////////////////////////////////////////////
// Base class for binary operators performing: Node * Node -> Node 
/////////////////////////////////////////////////////////////////////////////////////////

class base_NNtoN : public MEDDLY::binary_operation {
public:
    base_NNtoN(MEDDLY::binary_opname* opcode, MEDDLY::forest* arg1,
               MEDDLY::forest* arg2, MEDDLY::forest* res);

    virtual void computeDDEdge(const MEDDLY::dd_edge &ar1, const MEDDLY::dd_edge &ar2, 
                               MEDDLY::dd_edge &res, bool userFlag) override;

protected:
    virtual MEDDLY::node_handle compute(MEDDLY::node_handle a, MEDDLY::node_handle b) = 0;
    
    inline MEDDLY::ct_entry_key* 
    findResult(MEDDLY::node_handle a, MEDDLY::node_handle b, MEDDLY::node_handle &c);

    inline void 
    saveResult(MEDDLY::ct_entry_key* key, 
               MEDDLY::node_handle a, MEDDLY::node_handle b, MEDDLY::node_handle c);

    // utility to perform node unions
    MEDDLY::binary_operation* mddUnion;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Base class for unary operators performing: Node -> Node 
/////////////////////////////////////////////////////////////////////////////////////////

class base_NtoN : public MEDDLY::unary_operation {
public:
    base_NtoN(MEDDLY::unary_opname* opcode, MEDDLY::forest* arg1,
              MEDDLY::forest* res);

    virtual void computeDDEdge(const MEDDLY::dd_edge &ar1, 
                               MEDDLY::dd_edge &res, bool userFlag) override;

protected:
    virtual MEDDLY::node_handle compute(MEDDLY::node_handle a) = 0;
    
    inline MEDDLY::ct_entry_key* 
    findResult(MEDDLY::node_handle a, MEDDLY::node_handle &c);

    inline void 
    saveResult(MEDDLY::ct_entry_key* key, 
               MEDDLY::node_handle a, MEDDLY::node_handle c);

    // utility to perform node unions
    MEDDLY::binary_operation* mddUnion;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Base class for unary operators performing: Node -> int 
/////////////////////////////////////////////////////////////////////////////////////////

class base_NtoI : public MEDDLY::unary_operation {
public:
    base_NtoI(MEDDLY::unary_opname* opcode, MEDDLY::forest* arg1);

    virtual void computeDDEdge(const MEDDLY::dd_edge &ar1, 
                               MEDDLY::dd_edge &res, bool userFlag) override;

    void computeDDEdge(const MEDDLY::dd_edge &ar1, 
                       int &res, bool userFlag=false);

protected:
    virtual int compute(MEDDLY::node_handle a) = 0;
    
    inline MEDDLY::ct_entry_key* 
    findResult(MEDDLY::node_handle a, int &c);

    inline void 
    saveResult(MEDDLY::ct_entry_key* key, MEDDLY::node_handle a, int c);
};

/////////////////////////////////////////////////////////////////////////////////////////
// Base of unary operations with signature: node * integer -> node
/////////////////////////////////////////////////////////////////////////////////////////

class base_NItoN : public MEDDLY::operation {
public:
    base_NItoN(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
               MEDDLY::forest* _resF);
    virtual ~base_NItoN();    

    void computeDDEdge(const MEDDLY::dd_edge &ar, const int b, MEDDLY::dd_edge &res);

    bool checkForestCompatibility() const override;

protected:
    MEDDLY::forest* argF;
    MEDDLY::forest* resF;
    // opnd_type resultType;

    virtual MEDDLY::node_handle compute(MEDDLY::node_handle a, const int b) = 0;
    
    inline MEDDLY::ct_entry_key* 
    findResult(MEDDLY::node_handle a, const int b, MEDDLY::node_handle &c);

    inline void saveResult(MEDDLY::ct_entry_key* key, 
                           MEDDLY::node_handle a, const int b, MEDDLY::node_handle c);

    MEDDLY::binary_operation* mddUnion;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Base class for binary operators performing: Node * Node * Integer -> Node 
/////////////////////////////////////////////////////////////////////////////////////////

class base_NNItoN : public MEDDLY::operation {
public:
    base_NNItoN(MEDDLY::opname* opcode, MEDDLY::forest* arg1,
                MEDDLY::forest* arg2, MEDDLY::forest* res);
    virtual ~base_NNItoN();

    bool checkForestCompatibility() const override;

    // virtual void computeDDEdge(const MEDDLY::dd_edge &ar1, const MEDDLY::dd_edge &ar2, 
    //                            MEDDLY::dd_edge &res, bool userFlag) override;

protected:
    MEDDLY::forest* arg1F, *arg2F;
    MEDDLY::forest* resF;

    inline MEDDLY::ct_entry_key* 
    findResult(MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
               MEDDLY::node_handle &c);

    inline void 
    saveResult(MEDDLY::ct_entry_key* key, 
               //MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
               MEDDLY::node_handle c);

    // utility to perform node unions
    MEDDLY::binary_operation* mddUnion;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Base of unary operations with signature: node * integer -> node * node
/////////////////////////////////////////////////////////////////////////////////////////

class base_NItoNN : public MEDDLY::operation {
public:
    base_NItoNN(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
                MEDDLY::forest* _res1F, MEDDLY::forest* _res2F);
    virtual ~base_NItoNN();    

    void computeDDEdge(const MEDDLY::dd_edge &ar, const int b, 
                       MEDDLY::dd_edge &res1, MEDDLY::dd_edge &res2);

    bool checkForestCompatibility() const override;

protected:
    MEDDLY::forest* argF;
    MEDDLY::forest* res1F;
    MEDDLY::forest* res2F;

    virtual std::pair<MEDDLY::node_handle, MEDDLY::node_handle>
    compute(MEDDLY::node_handle a, const int b) = 0;
    
    inline MEDDLY::ct_entry_key* 
    findResult(MEDDLY::node_handle a, const int b, 
               std::pair<MEDDLY::node_handle, MEDDLY::node_handle> &c);

    inline void saveResult(MEDDLY::ct_entry_key* key, 
                           //MEDDLY::node_handle a, const int b, 
                           std::pair<MEDDLY::node_handle, MEDDLY::node_handle> c);

    MEDDLY::binary_operation* mddUnion;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Base class for binary operators performing: Node * Node * Integer -> Node * Node
/////////////////////////////////////////////////////////////////////////////////////////

class base_NNItoNN : public MEDDLY::operation {
public:
    base_NNItoNN(MEDDLY::opname* opcode, 
                 MEDDLY::forest* arg1, MEDDLY::forest* arg2, 
                 MEDDLY::forest* res1, MEDDLY::forest* res2);
    virtual ~base_NNItoNN();

    bool checkForestCompatibility() const override;

protected:
    MEDDLY::forest* arg1F, *arg2F;
    MEDDLY::forest* res1F, *res2F;

    virtual std::pair<MEDDLY::node_handle, MEDDLY::node_handle>
    compute(MEDDLY::node_handle a, MEDDLY::node_handle b, const int i) = 0;

    inline MEDDLY::ct_entry_key* 
    findResult(MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
               std::pair<MEDDLY::node_handle, MEDDLY::node_handle> &c);

    inline void 
    saveResult(MEDDLY::ct_entry_key* key, 
               //MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
               std::pair<MEDDLY::node_handle, MEDDLY::node_handle> c);

    // utility to perform node unions
    MEDDLY::binary_operation* mddUnion;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Base class for binary operators performing: Node * Node * Integer -> Node * Node * Node
/////////////////////////////////////////////////////////////////////////////////////////

class base_NNItoNNN : public MEDDLY::operation {
public:
    base_NNItoNNN(MEDDLY::opname* opcode, 
                  MEDDLY::forest* arg1, MEDDLY::forest* arg2, 
                  MEDDLY::forest* res1, MEDDLY::forest* res2, 
                  MEDDLY::forest* res3);
    virtual ~base_NNItoNNN();

    bool checkForestCompatibility() const override;

protected:
    MEDDLY::forest* arg1F, *arg2F;
    MEDDLY::forest* res1F, *res2F, *res3F;

    virtual std::tuple<MEDDLY::node_handle, MEDDLY::node_handle, MEDDLY::node_handle>
    compute(MEDDLY::node_handle a, MEDDLY::node_handle b, const int i) = 0;

    inline MEDDLY::ct_entry_key* 
    findResult(MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
               std::tuple<MEDDLY::node_handle, MEDDLY::node_handle, MEDDLY::node_handle> &c);

    inline void 
    saveResult(MEDDLY::ct_entry_key* key, 
               //MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
               std::tuple<MEDDLY::node_handle, MEDDLY::node_handle, MEDDLY::node_handle> c);
};

/////////////////////////////////////////////////////////////////////////////////////////













/////////////////////////////////////////////////////////////////////////////////////////
// S-Vectors
/////////////////////////////////////////////////////////////////////////////////////////

// TODO: si può eliminare NEG e avere solo decided/undecided, per poi scendere
// ricorsivamente con (p,q) o (q,p) per decidere il segno.
// NOTA2: la somma (A+B) può essere commutativa.
enum sv_sign { SVS_UNDECIDED, SVS_POS, SVS_NEG, SVS_TOTAL };

enum class ab_sum_t { A_PLUS_B, A_MINUS_B };

enum cmp_sign { CMP_UNDECIDED, CMP_POS, CMP_NEG, CMP_TOTAL };

/////////////////////////////////////////////////////////////////////////////////////////

class s_vectors : public base_NNItoN {
public:
    s_vectors(MEDDLY::opname* opcode, MEDDLY::forest* arg1,
              MEDDLY::forest* arg2, MEDDLY::forest* res,
              const variable_order *pivot_order);

    virtual MEDDLY::node_handle compute(MEDDLY::node_handle a, MEDDLY::node_handle b, 
                                        const bool is_potentially_conformant, 
                                        const ab_sum_t sum_or_diff, 
                                        const sv_sign sign_of_sum, const size_t lambda);

    void 
    computeDDEdge(const MEDDLY::dd_edge &ar1, const MEDDLY::dd_edge &ar2, 
                  const bool is_potentially_conformant, const ab_sum_t sum_or_diff, 
                  const sv_sign sign_of_sum, const size_t lambda,
                  MEDDLY::dd_edge &res);

protected:
    const variable_order *pivot_order; // pivoting order when proceeding by levels
};

// Factory of s_vectors operators for specific MDD forests
class s_vectors_opname : public MEDDLY::opname {
public:
    inline s_vectors_opname() 
    /**/ : MEDDLY::opname("S-Vectors") {}

    s_vectors* 
    buildOperation(MEDDLY::forest* arg1,
                   MEDDLY::forest* arg2, 
                   MEDDLY::forest* res,
                   const variable_order *pivot_order);
};

// /////////////////////////////////////////////////////////////////////////////////////////
// // Less Equal Squared Operator
// // Get the elements using the less-equal-but-not-equal-squared operator
// /////////////////////////////////////////////////////////////////////////////////////////

// class lesseq_sq : public base_NNtoN {
// public:
//     lesseq_sq(MEDDLY::binary_opname* opcode, MEDDLY::forest* arg1,
//               MEDDLY::forest* arg2, MEDDLY::forest* res,
//               const lesseq_sq_table* tab,
//               bool isEq, bool isB0,
//               bool _subtract, size_t lambda);

//     virtual MEDDLY::node_handle compute(MEDDLY::node_handle a, MEDDLY::node_handle b) override;

// protected:

//     const lesseq_sq_table* p_table; // operator's variations

//     const bool isPotentiallyEqual; // is a==b ?
//     const bool isBZero; // is b potentially zero?
//     const bool subtract; // return (a-b) instead of (a)
//     const size_t lambda; // is level-specific (>0) or not (=0)
// };

// // Factory of lesseq_sq operators for specific MDD forests
// class lesseq_sq_opname : public MEDDLY::binary_opname {
// public:
//     inline lesseq_sq_opname() 
//     /**/ : MEDDLY::binary_opname("LessEqSquared") {}

//     virtual MEDDLY::binary_operation* 
//     buildOperation(MEDDLY::forest* a1, 
//                    MEDDLY::forest* a2, 
//                    MEDDLY::forest* r) override;

//     lesseq_sq* 
//     buildOperation(MEDDLY::forest* a1, 
//                    MEDDLY::forest* a2, 
//                    MEDDLY::forest* r,
//                    const lesseq_sq_table* tab,
//                    bool isEq, bool isB0, bool subtract,
//                    size_t lambda);
// };

// // Table with all parametric s_vector instances
// class lesseq_sq_table {
//     // level * isEq * isB0 * subtract
//     mutable std::vector<std::vector<std::vector<std::vector<lesseq_sq*>>>> table;
// public:
//     lesseq_sq_table(MEDDLY::forest* forest, const variable_order *pivot_order);

//     lesseq_sq* get_op(size_t level, bool isPotentiallyEqual, 
//                       bool isBZero, bool subtract) const; 

//     MEDDLY::forest* forest;
//     const variable_order *pivot_order; // pivoting order when proceeding by levels
// };

/////////////////////////////////////////////////////////////////////////////////////////
// Less-Equal-but-Not-Equal-Squared Operator
// Get the elements using the less-equal-but-not-equal-squared operator
/////////////////////////////////////////////////////////////////////////////////////////

class leq_neq_sq : public base_NNItoN {
public:
    leq_neq_sq(MEDDLY::opname* opcode, MEDDLY::forest* forestMDD,
               const variable_order *pivot_order, const bool subtract);

    void 
    computeDDEdge(const MEDDLY::dd_edge &a, const MEDDLY::dd_edge &b, 
                  const bool is_potentially_equal, 
                  const bool is_b_potentially_zero,
                  const size_t lambda,
                  MEDDLY::dd_edge &res);
protected:

    MEDDLY::node_handle compute(MEDDLY::node_handle a, MEDDLY::node_handle b, int flags);

    const variable_order *pivot_order;
    const bool subtract; // compute (a-b) instead of (a)
};

// Factory of leq_neq_sq operators for specific MDD forests
class leq_neq_sq_opname : public MEDDLY::opname {
public:
    inline leq_neq_sq_opname() 
    /**/ : MEDDLY::opname("LessEqNeqSquared") {}

    leq_neq_sq* 
    buildOperation(MEDDLY::forest* forestMDD,
                   const variable_order *pivot_order,
                   const bool subtract);
};

/////////////////////////////////////////////////////////////////////////////////////////
// Reduction of elements that are less-equal-squared-but-not-equal up to lambda
/////////////////////////////////////////////////////////////////////////////////////////

class reduce : public base_NNItoNN {
public:
    reduce(MEDDLY::opname* opcode, MEDDLY::forest* forestMDD,
           const variable_order *pivot_order,
           const bool compute_differences);

    void 
    computeDDEdge(const MEDDLY::dd_edge &a, const MEDDLY::dd_edge &b, 
                  const bool is_potentially_equal, 
                  const bool is_b_potentially_zero,
                  const sv_sign sign_of_sum, 
                  const cmp_sign sign_of_comparison,
                  const size_t lambda,
                  MEDDLY::dd_edge &irreducibles, 
                  MEDDLY::dd_edge &reduced);

    // long counter_steps = 0; // TODO: remove

protected:
    virtual std::pair<MEDDLY::node_handle, MEDDLY::node_handle>
    compute(MEDDLY::node_handle a, MEDDLY::node_handle b, const int i) override;

    const variable_order *pivot_order; // pivoting order when proceeding by levels

    // utility to perform node operations
    MEDDLY::binary_operation *mddUnion;
    const bool compute_differences;
};

// Factory of reduce operators for specific MDD forests
class reduce_opname : public MEDDLY::opname {
public:
    inline reduce_opname() 
    /**/ : MEDDLY::opname("Reduce_LeqSq") {}

    reduce* 
    buildOperation(MEDDLY::forest *forestMDD,
                   const variable_order *pivot_order,
                   const bool compute_differences);
};

// /////////////////////////////////////////////////////////////////////////////////////////
// // Completion Procedure
// // Remove elements that are equal between [0..j-1] and smaller at var j
// /////////////////////////////////////////////////////////////////////////////////////////

// class compl_proc : public base_NNtoN {
// public:
//     compl_proc(MEDDLY::binary_opname* opcode, MEDDLY::forest* arg1,
//                MEDDLY::forest* arg2, MEDDLY::forest* res,
//                const compl_proc_table* tab,
//                size_t _rlvl);

//     virtual MEDDLY::node_handle compute(MEDDLY::node_handle a, MEDDLY::node_handle b) override;

// protected:

//     const compl_proc_table* p_table; // operator's variations

//     const size_t restricted_level; // completion level
// };

// // Factory of compl_proc operators for specific MDD forests
// class compl_proc_opname : public MEDDLY::binary_opname {
// public:
//     inline compl_proc_opname() 
//     /**/ : MEDDLY::binary_opname("CompletionProcedure") {}

//     virtual MEDDLY::binary_operation* 
//     buildOperation(MEDDLY::forest* a1, 
//                    MEDDLY::forest* a2, 
//                    MEDDLY::forest* r) override;

//     compl_proc* 
//     buildOperation(MEDDLY::forest* a1, 
//                    MEDDLY::forest* a2, 
//                    MEDDLY::forest* r,
//                    const compl_proc_table* tab,
//                    size_t restricted_level=0);
// };

// // Table with all parametric s_vector instances
// class compl_proc_table {
//     // level
//     std::vector<compl_proc*> table;
// public:
//     compl_proc_table(MEDDLY::forest* forest);

//     inline compl_proc* get_op(size_t level) const 
//     { return table[level]; }
// };

/////////////////////////////////////////////////////////////////////////////////////////
// Canonicalize the sign of a MDD
/////////////////////////////////////////////////////////////////////////////////////////

// Return the set where only one vector between each pair (v, -v) is kept
class sign_canon_mdd_op : public base_NtoN {
public:
    sign_canon_mdd_op(MEDDLY::unary_opname* oc, MEDDLY::forest* _argF,
                      MEDDLY::forest* _resF);
    virtual ~sign_canon_mdd_op();
public:
    virtual MEDDLY::node_handle compute(MEDDLY::node_handle a) override;
};

// Factory
class sign_canon_mdd_opname : public MEDDLY::unary_opname {
public:
    sign_canon_mdd_opname();
    
    virtual sign_canon_mdd_op* 
    buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF);
};

///////////////////////////////////////////////////////////////////////////////////////
// Look-up table of integer sets
///////////////////////////////////////////////////////////////////////////////////////

class LUT_int_set {
    std::vector< std::vector<int> > LUT;
    std::map< std::vector<int>, int > invLUT;
public:
    inline const std::vector<int>& 
    look_up(size_t i) const 
    { return LUT[i]; }

    inline size_t 
    insert(std::vector<int>&& key) {
        auto it = invLUT.find(key);
        if (it == invLUT.end()) { // new
            size_t new_pos = LUT.size();
            LUT.push_back(key);
            invLUT.insert(make_pair(std::move(key), new_pos));
            return new_pos;
        }
        else return it->second; // existing
    }
};

/////////////////////////////////////////////////////////////////////////////////////////
// Divisors of DD paths
/////////////////////////////////////////////////////////////////////////////////////////

// Returns the set of divisors of all descending paths of a node
class divisors_finder_mdd_op : public base_NtoI {
public:
    typedef int lut_key;
    divisors_finder_mdd_op(MEDDLY::unary_opname* oc, MEDDLY::forest* arg);

    inline const std::vector<int>& look_up(lut_key i) const { return lut.look_up(i); }
private:
    // Look-up table encoding the sets
    LUT_int_set lut;

protected:
    virtual lut_key compute(MEDDLY::node_handle a) override;
};

/////////////////////////////////////////////////////////////////////////////////////////

// Factory
class divisors_finder_mdd_opname : public MEDDLY::unary_opname {
public:
    divisors_finder_mdd_opname();
    virtual MEDDLY::unary_operation* buildOperation(MEDDLY::forest* ar);
};

/////////////////////////////////////////////////////////////////////////////////////////
// Divide vectors, when possible
// Divide all paths that are divisible by a given number
/////////////////////////////////////////////////////////////////////////////////////////

class vdivide_mdd_op : public base_NItoNN {
public:
    vdivide_mdd_op(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
                  MEDDLY::forest* _resF);
    virtual ~vdivide_mdd_op();    

protected:
    virtual std::pair<MEDDLY::node_handle, MEDDLY::node_handle>
    compute(MEDDLY::node_handle a, const int divisor) override;
};

/////////////////////////////////////////////////////////////////////////////////////////

// Factory
class vdivide_mdd_opname : public MEDDLY::opname {
public:
    vdivide_mdd_opname();
    virtual MEDDLY::operation* buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF);
};

/////////////////////////////////////////////////////////////////////////////////////////
// Multiply vectors by a constant
/////////////////////////////////////////////////////////////////////////////////////////

class vmult_op : public base_NItoN {
public:
    vmult_op(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
             MEDDLY::forest* _resF);
    virtual ~vmult_op();    

protected:
    virtual MEDDLY::node_handle 
    compute(MEDDLY::node_handle a, const int multiplier) override;

    friend class sign_canon_mdd_op;
};

/////////////////////////////////////////////////////////////////////////////////////////

// Factory
class vmult_opname : public MEDDLY::opname {
public:
    vmult_opname();
    virtual MEDDLY::operation* buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF);
};

/////////////////////////////////////////////////////////////////////////////////////////
// Included support
// Get the elements for which there is another one with smaller support
//  { a \in A | exists b \in B, b!=a, supp(b) included in supp(a) }
/////////////////////////////////////////////////////////////////////////////////////////

class support_inclusion_op : public base_NNtoN {
public:
    support_inclusion_op(MEDDLY::binary_opname* opcode, MEDDLY::forest* arg1,
                         MEDDLY::forest* arg2, MEDDLY::forest* res,
                         const support_inclusion_table* p_table, bool is_pot_eq_supp, 
                         bool subtract, size_t lambda);

    virtual MEDDLY::node_handle compute(MEDDLY::node_handle a, MEDDLY::node_handle b) override;

protected:
    const support_inclusion_table* p_table; // operator's variations

    const bool is_pot_eq_supp; // support is potentially equal
    const bool subtract; // return (a-b) instead of (a)
    const size_t lambda; // is level-specific (>0) or not (=0)
};

// Factory of support_inclusion_op operators for specific MDD forests
class support_inclusion_opname : public MEDDLY::binary_opname {
public:
    inline support_inclusion_opname() 
    /**/ : MEDDLY::binary_opname("SupportInclusion") {}

    virtual MEDDLY::binary_operation* 
    buildOperation(MEDDLY::forest* a1, 
                   MEDDLY::forest* a2, 
                   MEDDLY::forest* r) override;

    support_inclusion_op* 
    buildOperation(MEDDLY::forest* a1, 
                   MEDDLY::forest* a2, 
                   MEDDLY::forest* r,
                   const support_inclusion_table* p_table,
                   bool is_pot_eq_supp, bool subtract, size_t lambda);
};

// Table with all parametric s_vector instances
class support_inclusion_table {
    // level * isEq * subtract
    mutable std::vector<std::vector<std::vector<support_inclusion_op*>>> table;
    MEDDLY::forest* forest;
public:
    support_inclusion_table(MEDDLY::forest* forest, const variable_order *pivot_order);

    support_inclusion_op* get_op(size_t level, bool isPotentiallyEqual, bool subtract) const;

    const variable_order *pivot_order; // pivoting order when proceeding by levels
};

/////////////////////////////////////////////////////////////////////////////////////////
// Level/Value selector
// Selector for a subset of an MDD having a given value at a specified level lambda
/////////////////////////////////////////////////////////////////////////////////////////

class lv_selector_op : public base_NItoN {
public:
    lv_selector_op(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
                   MEDDLY::forest* _resF, const size_t lambda);
    virtual ~lv_selector_op();    

protected:
    virtual MEDDLY::node_handle 
    compute(MEDDLY::node_handle a, const int value) override;

    const size_t lambda; // level at which the operator performs the selection
};

/////////////////////////////////////////////////////////////////////////////////////////

// Factory
class lv_selector_opname : public MEDDLY::unary_opname {
public:
    lv_selector_opname();
    // virtual MEDDLY::operation* buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF) const;

    lv_selector_op* 
    buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF, const size_t lambda);
};

// Table with all parametric s_vector instances
class lv_selector_table {
    // level 
   mutable std::vector<lv_selector_op*> table;
   MEDDLY::forest* forest;
public:
    lv_selector_table(MEDDLY::forest* forest);

    lv_selector_op* get_op(size_t level) const;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Enumerate all the distinct values that appear at a specified MDD level
/////////////////////////////////////////////////////////////////////////////////////////

class domain_values_enumerator {
    size_t level;
    std::vector<bool> pos; // >= 0
    std::vector<bool> neg; // < 0
    std::vector<bool> visited;

    void visit(const MEDDLY::forest *argF, const MEDDLY::node_handle a);

public:
    domain_values_enumerator(size_t level) 
    /**/ : level(level), visited(1024), pos(4), neg(4) { }

    void visit(const MEDDLY::dd_edge& dd) { 
        visit(static_cast<const MEDDLY::forest*>(dd.getForest()), dd.getNode()); 
    }

    void get_values(std::vector<int>& out) const;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Determine the smallest degree in a set, by level
/////////////////////////////////////////////////////////////////////////////////////////

enum class degree_type {
    BY_VALUE, BY_SUPPORT
};
inline int get_degree_of(int v, degree_type degtype) {
    assert(v >= 0);
    return (degtype == degree_type::BY_VALUE) ? v : (v>0);
}

/////////////////////////////////////////////////////////////////////////////////////////

class smallest_degree_op : public base_NtoI {
public:
    smallest_degree_op(MEDDLY::unary_opname* oc, 
                       MEDDLY::forest* arg,
                       const smallest_degree_table* p_table,
                       const degree_type degtype,
                       const size_t lambda);

protected:
    const smallest_degree_table* p_table; // operator's variations
    const degree_type degtype;
    const size_t lambda; // is level-specific (>0) or not (=0)

    virtual int compute(MEDDLY::node_handle a) override;
};

/////////////////////////////////////////////////////////////////////////////////////////

// Factory
class smallest_degree_opname : public MEDDLY::unary_opname {
public:
    smallest_degree_opname();
    // virtual MEDDLY::unary_operation* buildOperation(MEDDLY::forest* ar) const override;
    smallest_degree_op* buildOperation(MEDDLY::forest* ar, const smallest_degree_table* p_table, 
                                       const degree_type degtype, size_t lambda);
};

// Table with all parametric s_vector instances
class smallest_degree_table {
    // level * degtype
    mutable std::vector<std::vector<smallest_degree_op*>> table;
    MEDDLY::forest* forest;
public:
    smallest_degree_table(MEDDLY::forest* forest, const variable_order *pivot_order);

    smallest_degree_op* get_op(size_t level, const degree_type degtype) const;

    const variable_order *pivot_order; // pivoting order when proceeding by levels
};

/////////////////////////////////////////////////////////////////////////////////////////
// Enumerate all degrees in a set, by level
/////////////////////////////////////////////////////////////////////////////////////////

class degree_finder_op : public base_NtoI {
public:
    degree_finder_op(MEDDLY::unary_opname* oc, 
                     MEDDLY::forest* arg,
                     const degree_finder_table* p_table,
                     const degree_type degtype,
                     const size_t lambda);

    typedef int lut_key;

protected:
    const degree_finder_table* p_table; // operator's variations
    const degree_type degtype;
    const size_t lambda; // is level-specific (>0) or not (=0)

    virtual int compute(MEDDLY::node_handle a) override;
};

/////////////////////////////////////////////////////////////////////////////////////////

// Factory
class degree_finder_opname : public MEDDLY::unary_opname {
public:
    degree_finder_opname();
    degree_finder_op* buildOperation(MEDDLY::forest* ar, const degree_finder_table* p_table, 
                                     const degree_type degtype, size_t lambda);
};

// Table with all parametric s_vector instances
class degree_finder_table {
    // level * degree_type
    mutable std::vector<std::vector<degree_finder_op*>> table;
    MEDDLY::forest* forest;
public:
    degree_finder_table(MEDDLY::forest* forest, const variable_order *pivot_order);

    degree_finder_op* get_op(size_t level, const degree_type degtype) const;

    const variable_order *pivot_order; // pivoting order when proceeding by levels

    inline const std::vector<int>& look_up(degree_finder_op::lut_key i) const 
    { return lut.look_up(i); }

private:
    // Look-up table encoding the sets
    mutable LUT_int_set lut;
    friend class degree_finder_op;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Select all elements that have a specific degree
/////////////////////////////////////////////////////////////////////////////////////////

class degree_selector_op : public base_NItoN {
public:
    degree_selector_op(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
                       MEDDLY::forest* _resF,
                       const degree_selector_table* p_table,
                       const degree_type degtype,
                       const size_t lambda);
    virtual ~degree_selector_op();    

protected:
    const degree_selector_table* p_table; // operator's variations
    const degree_type degtype;
    const size_t lambda; // is level-specific (>0) or not (=0)

protected:
    virtual MEDDLY::node_handle 
    compute(MEDDLY::node_handle a, const int degree) override;
};

/////////////////////////////////////////////////////////////////////////////////////////

// Factory
class degree_selector_opname : public MEDDLY::unary_opname {
public:
    degree_selector_opname();

    degree_selector_op* 
    buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF,
                   const degree_selector_table* p_table, 
                   const degree_type degtype, size_t lambda);
};

// Table with all parametric s_vector instances
class degree_selector_table {
    // level 
    mutable std::vector<std::vector<degree_selector_op*>> table;
    MEDDLY::forest* forest;
public:
    degree_selector_table(MEDDLY::forest* forest, const variable_order *pivot_order);

    degree_selector_op* get_op(size_t level, const degree_type degtype) const;

    const variable_order *pivot_order; // pivoting order when proceeding by levels
};

/////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////
#endif // __OPERATORS_H__
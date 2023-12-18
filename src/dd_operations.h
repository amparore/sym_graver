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

class s_vectors2_opname;
class s_vectors2;

DD_EXTERN s_vectors2_opname *S_VECTORS2_OPNAME;
DD_EXTERN s_vectors2 *S_VECTORS2;


class lesseq_sq_opname;
class lesseq_sq_table;
class lesseq_sq;

DD_EXTERN lesseq_sq_opname *LESSEQ_SQ_OPNAME;
DD_EXTERN lesseq_sq_table *LESSEQ_SQ_OPS;

// class compl_proc_opname;
// class compl_proc_table;
// class compl_proc;

// DD_EXTERN compl_proc_opname *COMPL_PROC_OPNAME;
// DD_EXTERN compl_proc_table *COMPL_PROC_OPS;

class sign_canon_mdd_opname;
class sign_canon_table;
DD_EXTERN sign_canon_mdd_opname* SIGN_CANON_OPNAME;
DD_EXTERN sign_canon_table* SIGN_CANON_OPS;

class vmult_opname;
class vmult_op;
DD_EXTERN vmult_opname* VMULT_OPNAME;
DD_EXTERN vmult_op* VMULT;


class vcanon_mdd_opname;
class divisors_finder_mdd_opname;

class vcanon_mdd_op;
class divisors_finder_mdd_op;

// DD_EXTERN lesseq_sq* LESSEQ_SQ;
// DD_EXTERN lesseq_sq* REDUCE_LESSEQ_SQ;

DD_EXTERN vcanon_mdd_opname* VCANON_DIVISORS_MDD_OPNAME;
DD_EXTERN vcanon_mdd_op* VCANON_DIVISORS_MDD;
DD_EXTERN vcanon_mdd_opname* VCANON_DIVIDE_MDD_OPNAME;
DD_EXTERN vcanon_mdd_op* VCANON_DIVIDE_MDD;

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
};
std::ostream& operator<< (std::ostream&, const mdd_printer&);

inline mdd_printer print_mdd(MEDDLY::dd_edge dd, const variable_order& vorder,
                             bool write_header=true, bool write_counts=false) {
    return mdd_printer{ .dd=dd, .vorder=vorder, .write_header=write_header, 
                        .write_counts=write_counts };
}

/////////////////////////////////////////////////////////////////////////////////////////

void check_level_bound(MEDDLY::forest* forest, int level, int bound);

MEDDLY::dd_edge 
mdd_from_vectors(const std::vector<std::vector<int>>& vecs, 
                 MEDDLY::forest* forestMDD, bool make_symmetric);

inline bool is_emptyset(const MEDDLY::dd_edge& e) {
    MEDDLY::expert_forest *forest = static_cast<MEDDLY::expert_forest *>(e.getForest());
    return e.getNode() == forest->handleForValue(false);
}

MEDDLY::dd_edge
selector_for_value_at_level(MEDDLY::forest *forestMDD, int value, size_t level);

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
                unsigned pos, MEDDLY::expert_forest* resF,  
                MEDDLY::binary_operation* mddUnion) ;

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
    base_NNtoN(MEDDLY::binary_opname* opcode, MEDDLY::expert_forest* arg1,
               MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res);

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
    base_NtoN(MEDDLY::unary_opname* opcode, MEDDLY::expert_forest* arg1,
              MEDDLY::expert_forest* res);

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
    base_NtoI(MEDDLY::unary_opname* opcode, MEDDLY::expert_forest* arg1);

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
    base_NItoN(MEDDLY::opname* opcode, MEDDLY::expert_forest* _argF,
               MEDDLY::expert_forest* _resF);
    virtual ~base_NItoN();    

    void computeDDEdge(const MEDDLY::dd_edge &ar, const int b, MEDDLY::dd_edge &res);

    bool checkForestCompatibility() const override;

protected:
    MEDDLY::expert_forest* argF;
    MEDDLY::expert_forest* resF;
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
    base_NNItoN(MEDDLY::opname* opcode, MEDDLY::expert_forest* arg1,
                MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res);

    bool checkForestCompatibility() const override;

    // virtual void computeDDEdge(const MEDDLY::dd_edge &ar1, const MEDDLY::dd_edge &ar2, 
    //                            MEDDLY::dd_edge &res, bool userFlag) override;

protected:
    MEDDLY::expert_forest* arg1F, *arg2F;
    MEDDLY::expert_forest* resF;

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













/////////////////////////////////////////////////////////////////////////////////////////
// S-Vectors
/////////////////////////////////////////////////////////////////////////////////////////

enum sv_sign { SVS_UNDECIDED, SVS_POS, SVS_NEG, SVS_TOTAL };

// class s_vectors : public base_NNtoN {
// public:
//     s_vectors(MEDDLY::binary_opname* opcode, MEDDLY::expert_forest* arg1,
//               MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res,
//               const s_vectors_table *tab, bool isConf, bool invertB, 
//               const sv_sign sign_of_sum, size_t lambda);

//     virtual MEDDLY::node_handle compute(MEDDLY::node_handle a, MEDDLY::node_handle b) override;
// protected:
//     const s_vectors_table* p_table; // s_vector variations
//     const bool isPotentiallyConformant; // summed vector signs are potentially conformant?
//     const bool invertB; // sum with b or -b
//     const sv_sign sign_of_sum; // should we invert the sum? is this decided yet?
//     const size_t lambda; // is level-specific (>0) or not (=0)
// };

// // Factory of s_vectors operators for specific MDD forests
// class s_vectors_opname : public MEDDLY::binary_opname {
// public:
//     inline s_vectors_opname() 
//     /**/ : MEDDLY::binary_opname("S-Vectors") {}

//     virtual MEDDLY::binary_operation* 
//     buildOperation(MEDDLY::expert_forest* arg1,
//                    MEDDLY::expert_forest* arg2, 
//                    MEDDLY::expert_forest* res) override;

//     s_vectors* 
//     buildOperation(MEDDLY::expert_forest* a1, 
//                    MEDDLY::expert_forest* a2, 
//                    MEDDLY::expert_forest* r,
//                    const s_vectors_table *tab, 
//                    bool isConf, bool invertB, 
//                    const sv_sign sign_of_sum, size_t lambda);
// };

// // Table with all parametric s_vector instances
// class s_vectors_table {
//     // level * isConf * invertB * sign_of_sum
//     mutable std::vector<std::vector<std::vector<std::vector<s_vectors*>>>> table; 
// public:
//     s_vectors_table(MEDDLY::expert_forest* forest, const variable_order *pivot_order);

//     s_vectors* get_op(size_t level, bool isPotentiallyConformant, 
//                       bool invertB, const sv_sign sign_of_sum) const;

//     MEDDLY::expert_forest* forest;
//     const variable_order *pivot_order; // pivoting order when proceeding by levels
// };

/////////////////////////////////////////////////////////////////////////////////////////
// S-Vectors ternary
/////////////////////////////////////////////////////////////////////////////////////////

class s_vectors2 : public base_NNItoN {
public:
    s_vectors2(MEDDLY::opname* opcode, MEDDLY::expert_forest* arg1,
               MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res,
               const variable_order *pivot_order);

    virtual MEDDLY::node_handle compute(MEDDLY::node_handle a, MEDDLY::node_handle b, 
                                        const bool is_potentially_conformant, 
                                        const bool invertB, 
                                        const sv_sign sign_of_sum, const size_t lambda);

    void 
    computeDDEdge(const MEDDLY::dd_edge &ar1, const MEDDLY::dd_edge &ar2, 
                  const bool is_potentially_conformant, const bool invertB, 
                  const sv_sign sign_of_sum, const size_t lambda,
                  MEDDLY::dd_edge &res);

protected:
    const variable_order *pivot_order; // pivoting order when proceeding by levels
};

// Factory of s_vectors2 operators for specific MDD forests
class s_vectors2_opname : public MEDDLY::opname {
public:
    inline s_vectors2_opname() 
    /**/ : MEDDLY::opname("S-Vectors") {}

    s_vectors2* 
    buildOperation(MEDDLY::expert_forest* arg1,
                   MEDDLY::expert_forest* arg2, 
                   MEDDLY::expert_forest* res,
                   const variable_order *pivot_order);
};

/////////////////////////////////////////////////////////////////////////////////////////
// Less Equal Squared Operator
// Get the elements using the less-equal-but-not-equal-squared operator
/////////////////////////////////////////////////////////////////////////////////////////

class lesseq_sq : public base_NNtoN {
public:
    lesseq_sq(MEDDLY::binary_opname* opcode, MEDDLY::expert_forest* arg1,
              MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res,
              const lesseq_sq_table* tab,
              bool isEq, bool isB0,
              bool _subtract, size_t lambda);

    virtual MEDDLY::node_handle compute(MEDDLY::node_handle a, MEDDLY::node_handle b) override;

protected:

    const lesseq_sq_table* p_table; // operator's variations

    const bool isPotentiallyEqual; // is a==b ?
    const bool isBZero; // is b potentially zero?
    const bool subtract; // return (a-b) instead of (a)
    const size_t lambda; // is level-specific (>0) or not (=0)
};

// Factory of lesseq_sq operators for specific MDD forests
class lesseq_sq_opname : public MEDDLY::binary_opname {
public:
    inline lesseq_sq_opname() 
    /**/ : MEDDLY::binary_opname("LessEqSquared") {}

    virtual MEDDLY::binary_operation* 
    buildOperation(MEDDLY::expert_forest* a1, 
                   MEDDLY::expert_forest* a2, 
                   MEDDLY::expert_forest* r) override;

    lesseq_sq* 
    buildOperation(MEDDLY::expert_forest* a1, 
                   MEDDLY::expert_forest* a2, 
                   MEDDLY::expert_forest* r,
                   const lesseq_sq_table* tab,
                   bool isEq, bool isB0, bool subtract,
                   size_t lambda);
};

// Table with all parametric s_vector instances
class lesseq_sq_table {
    // level * isEq * isB0 * subtract
    mutable std::vector<std::vector<std::vector<std::vector<lesseq_sq*>>>> table;
public:
    lesseq_sq_table(MEDDLY::expert_forest* forest, const variable_order *pivot_order);

    lesseq_sq* get_op(size_t level, bool isPotentiallyEqual, 
                      bool isBZero, bool subtract) const; 

    MEDDLY::expert_forest* forest;
    const variable_order *pivot_order; // pivoting order when proceeding by levels
};

// /////////////////////////////////////////////////////////////////////////////////////////
// // Completion Procedure
// // Remove elements that are equal between [0..j-1] and smaller at var j
// /////////////////////////////////////////////////////////////////////////////////////////

// class compl_proc : public base_NNtoN {
// public:
//     compl_proc(MEDDLY::binary_opname* opcode, MEDDLY::expert_forest* arg1,
//                MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res,
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
//     buildOperation(MEDDLY::expert_forest* a1, 
//                    MEDDLY::expert_forest* a2, 
//                    MEDDLY::expert_forest* r) override;

//     compl_proc* 
//     buildOperation(MEDDLY::expert_forest* a1, 
//                    MEDDLY::expert_forest* a2, 
//                    MEDDLY::expert_forest* r,
//                    const compl_proc_table* tab,
//                    size_t restricted_level=0);
// };

// // Table with all parametric s_vector instances
// class compl_proc_table {
//     // level
//     std::vector<compl_proc*> table;
// public:
//     compl_proc_table(MEDDLY::expert_forest* forest);

//     inline compl_proc* get_op(size_t level) const 
//     { return table[level]; }
// };

/////////////////////////////////////////////////////////////////////////////////////////
// Canonicalize the sign of a MDD
/////////////////////////////////////////////////////////////////////////////////////////

// Return the set where only one vector between each pair (v, -v) is kept
class sign_canon_mdd_op : public base_NtoN {
public:
    sign_canon_mdd_op(MEDDLY::unary_opname* oc, MEDDLY::expert_forest* _argF,
                      MEDDLY::expert_forest* _resF, const sign_canon_table* p_table,
                      bool isZero);
    virtual ~sign_canon_mdd_op();
public:
    virtual MEDDLY::node_handle compute(MEDDLY::node_handle a) override;

protected:
    const sign_canon_table* p_table; // operator's variations
    const bool isZero; // is following a zero path (true) or is it inverting a vector (false)?
};

/////////////////////////////////////////////////////////////////////////////////////////

// Factory
class sign_canon_mdd_opname : public MEDDLY::unary_opname {
public:
    sign_canon_mdd_opname();
    
    virtual sign_canon_mdd_op* 
    buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF, 
                   const sign_canon_table* p_table, bool isZero);
};
// Table with all parametric sign_canon_mdd_op instances
class sign_canon_table {
    // isZero
    mutable std::vector<sign_canon_mdd_op*> table;
    MEDDLY::expert_forest* forest;
public:
    sign_canon_table(MEDDLY::expert_forest* forest);

    sign_canon_mdd_op* get_op(bool isZero) const;
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
    divisors_finder_mdd_op(MEDDLY::unary_opname* oc, MEDDLY::expert_forest* arg);

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
// Canonicalize vectors
// Divide all paths that are divisible by a given number
/////////////////////////////////////////////////////////////////////////////////////////

class vcanon_mdd_op : public base_NItoN {
public:
    vcanon_mdd_op(MEDDLY::opname* opcode, MEDDLY::expert_forest* _argF,
                  MEDDLY::expert_forest* _resF, bool _divide);
    virtual ~vcanon_mdd_op();    

protected:
    virtual MEDDLY::node_handle 
    compute(MEDDLY::node_handle a, const int divisor) override;

    // should divide the elements by the divisor?
    bool divide;
};

/////////////////////////////////////////////////////////////////////////////////////////

// Factory
class vcanon_mdd_opname : public MEDDLY::unary_opname {
public:
    vcanon_mdd_opname(bool _divide);
    virtual MEDDLY::operation* buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF);
private:
    bool divide;
};

/////////////////////////////////////////////////////////////////////////////////////////
// Multiply vectors by a constant
/////////////////////////////////////////////////////////////////////////////////////////

class vmult_op : public base_NItoN {
public:
    vmult_op(MEDDLY::opname* opcode, MEDDLY::expert_forest* _argF,
             MEDDLY::expert_forest* _resF);
    virtual ~vmult_op();    

protected:
    virtual MEDDLY::node_handle 
    compute(MEDDLY::node_handle a, const int multiplier) override;
};

/////////////////////////////////////////////////////////////////////////////////////////

// Factory
class vmult_opname : public MEDDLY::unary_opname {
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
    support_inclusion_op(MEDDLY::binary_opname* opcode, MEDDLY::expert_forest* arg1,
                         MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res,
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
    buildOperation(MEDDLY::expert_forest* a1, 
                   MEDDLY::expert_forest* a2, 
                   MEDDLY::expert_forest* r) override;

    support_inclusion_op* 
    buildOperation(MEDDLY::expert_forest* a1, 
                   MEDDLY::expert_forest* a2, 
                   MEDDLY::expert_forest* r,
                   const support_inclusion_table* p_table,
                   bool is_pot_eq_supp, bool subtract, size_t lambda);
};

// Table with all parametric s_vector instances
class support_inclusion_table {
    // level * isEq * subtract
    mutable std::vector<std::vector<std::vector<support_inclusion_op*>>> table;
    MEDDLY::expert_forest* forest;
public:
    support_inclusion_table(MEDDLY::expert_forest* forest, const variable_order *pivot_order);

    support_inclusion_op* get_op(size_t level, bool isPotentiallyEqual, bool subtract) const;

    const variable_order *pivot_order; // pivoting order when proceeding by levels
};

/////////////////////////////////////////////////////////////////////////////////////////
// Level/Value selector
// Selector for a subset of an MDD having a given value at a specified level lambda
/////////////////////////////////////////////////////////////////////////////////////////

class lv_selector_op : public base_NItoN {
public:
    lv_selector_op(MEDDLY::opname* opcode, MEDDLY::expert_forest* _argF,
                   MEDDLY::expert_forest* _resF, const size_t lambda);
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
   MEDDLY::expert_forest* forest;
public:
    lv_selector_table(MEDDLY::expert_forest* forest);

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

    void visit(const MEDDLY::expert_forest *argF, const MEDDLY::node_handle a);

public:
    domain_values_enumerator(size_t level) 
    /**/ : level(level), visited(1024), pos(4), neg(4) { }

    void visit(const MEDDLY::dd_edge& dd) { 
        visit(static_cast<const MEDDLY::expert_forest*>(dd.getForest()), dd.getNode()); 
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
                       MEDDLY::expert_forest* arg,
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
    MEDDLY::expert_forest* forest;
public:
    smallest_degree_table(MEDDLY::expert_forest* forest, const variable_order *pivot_order);

    smallest_degree_op* get_op(size_t level, const degree_type degtype) const;

    const variable_order *pivot_order; // pivoting order when proceeding by levels
};

/////////////////////////////////////////////////////////////////////////////////////////
// Enumerate all degrees in a set, by level
/////////////////////////////////////////////////////////////////////////////////////////

class degree_finder_op : public base_NtoI {
public:
    degree_finder_op(MEDDLY::unary_opname* oc, 
                     MEDDLY::expert_forest* arg,
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
    MEDDLY::expert_forest* forest;
public:
    degree_finder_table(MEDDLY::expert_forest* forest, const variable_order *pivot_order);

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
    degree_selector_op(MEDDLY::opname* opcode, MEDDLY::expert_forest* _argF,
                       MEDDLY::expert_forest* _resF,
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
    MEDDLY::expert_forest* forest;
public:
    degree_selector_table(MEDDLY::expert_forest* forest, const variable_order *pivot_order);

    degree_selector_op* get_op(size_t level, const degree_type degtype) const;

    const variable_order *pivot_order; // pivoting order when proceeding by levels
};

/////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////
#endif // __OPERATORS_H__
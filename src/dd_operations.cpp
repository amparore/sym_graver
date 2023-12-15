#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cassert>

#define DD_OPERATIONS_IMPL 1
#include "dd_operations.h"

#include "math_utils.h"
#include "variable_order.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////

void account_col_space(std::vector<size_t>& col_spaces, const int *const data);

std::ostream& operator<< (std::ostream& os, const mdd_printer& mp)
{
    const size_t m = mp.dd.getForest()->getDomain()->getNumVariables();
    int values[m];

    std::vector<size_t> col_spaces(m);
    size_t num_rows = 0;
    for (MEDDLY::enumerator path(mp.dd); path != 0; ++path) {
        const int *assignments = path.getAssignments() + 1;
        for(size_t var=0; var < m; var++) {
            const size_t lvl = mp.vorder.var2lvl(var);
            values[var] = NodeToZ(assignments[lvl]);
        }
        account_col_space(col_spaces, values);
        ++num_rows;
    }

    // print counters
    if (mp.write_counts)
        os << num_rows << " " << m << endl;

    // print header
    if (mp.write_header) {
        for(size_t var=0; var < m; var++) {
            const size_t lvl = mp.vorder.var2lvl(var);
            const char* var_name = mp.dd.getForest()->getDomain()->getVar(lvl + 1)->getName();
            col_spaces[var] = max(col_spaces[var], strlen(var_name));
            os << setw(col_spaces[var]) << var_name << " ";
        }
        os << endl;
    }
    // print paths/tuples/rows
    for (MEDDLY::enumerator path(mp.dd); path != 0; ++path) {
        const int *assignments = path.getAssignments() + 1;
        for(size_t var=0; var < m; var++) {
            const size_t lvl = mp.vorder.var2lvl(var);
            os << setw(col_spaces[var]) << right << NodeToZ(assignments[lvl]) << " ";
        }
        os << endl;
    }
    return os;
}

/////////////////////////////////////////////////////////////////////////////////////////

void check_level_bound(MEDDLY::forest* forest, int level, int bound)
{
    // if (bound > 20)
    //     throw std::exception();

    MEDDLY::expert_domain* dom_exp = (MEDDLY::expert_domain*)forest->getDomain();
    if (dom_exp->getVariableBound(level) < bound) {
        // cout << "resizing bound of level "<<level<<" from "<<dom_exp->getVariableBound(level)
        //      << " to "<<bound<<endl;
        dom_exp->enlargeVariableBound(level, false, bound);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge 
mdd_from_vectors(const std::vector<std::vector<int>>& vecs, 
                 MEDDLY::forest* forestMDD, bool make_symmetric) 
{
    MEDDLY::expert_domain* dom_exp = (MEDDLY::expert_domain*)forestMDD->getDomain();
    const size_t num_rows = make_symmetric ? 2*vecs.size() : vecs.size();
    MEDDLY::dd_edge set(forestMDD);
    if (num_rows) {
        std::vector<std::vector<int>> ins_data(num_rows);
        std::vector<const int*> ins_ptr;
        ins_ptr.resize(num_rows);

        for (size_t i=0; i<num_rows; i++) {
            size_t row = i % vecs.size();
            int sign = (i < vecs.size()) ? 1 : -1;
            assert(dom_exp->getNumVariables() == vecs.at(row).size());
            ins_data[i].resize(dom_exp->getNumVariables() + 1);

            // encode integers and verify bounds
            for (size_t j=0; j<vecs[row].size(); j++) {
                ins_data[i][j+1] = ZtoNode(sign * vecs[row][j]); // encode integers
                if (ins_data[i][j+1] >= dom_exp->getVariableBound(j+1)) {
                    check_level_bound(forestMDD, j+1, ins_data[i][j+1] + 1);
                }
            }
            ins_ptr[i] = ins_data[i].data();
        }
        forestMDD->createEdge(ins_ptr.data(), ins_ptr.size(), set);
    }
    return set;
}


/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
selector_for_sign_at_level(MEDDLY::forest *forestMDD, size_t level, int sign, int threshold) 
{
    MEDDLY::expert_domain* dom_exp = (MEDDLY::expert_domain*)forestMDD->getDomain();

    int bound = dom_exp->getVariableBound(level);
    bool terms[bound+1];
    for (int i=0; i<bound+1; i++)
        terms[i] = (NodeToZ(i) * sign) >= threshold;

    MEDDLY::dd_edge level_sel(forestMDD);
    forestMDD->createEdgeForVar(level, false, terms, level_sel);
    return level_sel;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
selector_for_value_at_level(MEDDLY::forest *forestMDD, int value, size_t level) 
{
    MEDDLY::expert_domain* dom_exp = (MEDDLY::expert_domain*)forestMDD->getDomain();

    int bound = dom_exp->getVariableBound(level);
    bool terms[bound+1];
    for (int i=0; i<bound+1; i++)
        terms[i] = (NodeToZ(i) == value);

    MEDDLY::dd_edge level_sel(forestMDD);
    forestMDD->createEdgeForVar(level, false, terms, level_sel);
    return level_sel;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
selector_for_nonnegatives(MEDDLY::forest *forestMDD) 
{
    MEDDLY::dd_edge sel(forestMDD);
    for (int lvl=1; lvl<=forestMDD->getDomain()->getNumVariables(); lvl++) {
        MEDDLY::dd_edge level_sel = selector_for_sign_at_level(forestMDD, lvl, +1);
        if (lvl == 1)
            sel = level_sel;
        else
            MEDDLY::apply(MEDDLY::INTERSECTION, sel, level_sel, sel);
    }
    return sel;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
selector_for_nonnegatives(MEDDLY::forest *forestMDD, 
                          const std::vector<size_t> *rmnnz_levels, size_t rmnnz_step) 
{
    MEDDLY::dd_edge sel(forestMDD);
    for (int lvl=1; lvl<=forestMDD->getDomain()->getNumVariables(); lvl++) {
        if (rmnnz_levels==nullptr || (*rmnnz_levels)[lvl - 1] >= rmnnz_step) {
            // cout << "removing nonzeroes for level "<<lvl<<endl;
            MEDDLY::dd_edge level_sel = selector_for_sign_at_level(forestMDD, lvl, +1);
            if (is_emptyset(sel))
                sel = level_sel;
            else
                MEDDLY::apply(MEDDLY::INTERSECTION, sel, level_sel, sel);
        }
    }
    return sel;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
selector_for_zeros_up_to_level(MEDDLY::forest *forestMDD, size_t max_level)
{
    MEDDLY::dd_edge sel(forestMDD);
    for (int lvl=1; lvl<=max_level; lvl++) {
        MEDDLY::dd_edge level_sel = selector_for_value_at_level(forestMDD, 0, lvl);
        if (lvl == 1)
            sel = level_sel;
        else
            MEDDLY::apply(MEDDLY::INTERSECTION, sel, level_sel, sel);
    }
    return sel;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
zero_vector(MEDDLY::forest *forestMDD) {
    MEDDLY::expert_domain* dom_exp = (MEDDLY::expert_domain*)forestMDD->getDomain();
    MEDDLY::dd_edge zero(forestMDD);
    int num_levels = dom_exp->getNumVariables();
    std::vector<int> minterm(num_levels + 1);
    std::fill(minterm.begin(), minterm.end(), 0);
    const int* ins_ptr[1];
    ins_ptr[0] = minterm.data();
    forestMDD->createEdge(ins_ptr, 1, zero);
    return zero;
}

/////////////////////////////////////////////////////////////////////////////////////////

// Get the effective size (in node domain) of an unpacked_node
inline size_t 
get_node_size(MEDDLY::unpacked_node* n) {
    if (n->isFull()) {
        // find exact sizes (exclude trailing zeroes)
        size_t size = n->getSize();
        while (size > 0 && n->d(size - 1) == 0) 
            --size;
        return size;
    }
    else {
        // return the size containing the last nonzero
        return n->i( n->getNNZs()-1 ) + 1;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

void unionNodes(MEDDLY::unpacked_node* C, MEDDLY::node_handle node, 
                unsigned pos, MEDDLY::expert_forest* resF,  
                MEDDLY::binary_operation* mddUnion) 
{
    if (node == 0)
        return;
    assert(pos < C->getSize());
    assert(C->isFull());
    if (C->d(pos) == 0)
        C->d_ref(pos) = node; // resF->linkNode(node); // FIXME: check
    else { // perform union
        MEDDLY::dd_edge dd1(resF), dd2(resF), union_dd(resF);
        dd1.set(node);
        dd2.set(C->d(pos));
        mddUnion->computeTemp(dd1, dd2, union_dd);
        C->set_d(pos, union_dd);
        assert(resF->getNodeLevel(C->d(pos)) <= C->getLevel());
    }
}

/////////////////////////////////////////////////////////////////////////////////////////



















/////////////////////////////////////////////////////////////////////////////////////////
// Base class for binary operators performing: Node * Node -> Node 
/////////////////////////////////////////////////////////////////////////////////////////

base_NNtoN::base_NNtoN(MEDDLY::binary_opname* opcode, MEDDLY::expert_forest* arg1, 
                       MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res)
  : MEDDLY::binary_operation(opcode, 1, arg1, arg2, res)
{
    MEDDLY::ct_entry_type* et;
    et = new MEDDLY::ct_entry_type(opcode->getName(), "NN:N");
    et->setForestForSlot(0, arg1);
    et->setForestForSlot(1, arg2);
    et->setForestForSlot(3, res);
    registerEntryType(0, et);
    buildCTs();

    mddUnion = MEDDLY::getOperation(MEDDLY::UNION, res, res, res);
}

/////////////////////////////////////////////////////////////////////////////////////////
 
inline MEDDLY::ct_entry_key* 
base_NNtoN::findResult(MEDDLY::node_handle a, MEDDLY::node_handle b, MEDDLY::node_handle &c) 
{
    MEDDLY::ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    assert(CTsrch);
    if (can_commute && a > b) {
        CTsrch->writeN(b);
        CTsrch->writeN(a);
    } else {
        CTsrch->writeN(a);
        CTsrch->writeN(b);
    }
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    c = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
}

inline void 
base_NNtoN::saveResult(MEDDLY::ct_entry_key* key, 
                       MEDDLY::node_handle a, MEDDLY::node_handle b, MEDDLY::node_handle c) 
{
    CTresult[0].reset();
    CTresult[0].writeN(c);
    CT0->addEntry(key, CTresult[0]);
}

/////////////////////////////////////////////////////////////////////////////////////////

void base_NNtoN::computeDDEdge(const MEDDLY::dd_edge &ar1, const MEDDLY::dd_edge &ar2, 
                                       MEDDLY::dd_edge &res, bool userFlag)
{
    MEDDLY::node_handle cnode = compute(ar1.getNode(), ar2.getNode());
    // const int num_levels = resF->getDomain()->getNumVariables();
    // if (userFlag && resF->isQuasiReduced() && cnode != resF->getTransparentNode()
    //     && resF->getNodeLevel(cnode) < num_levels) 
    // {
    //     MEDDLY::node_handle temp = (resF)->makeNodeAtLevel(num_levels, cnode);
    //     resF->unlinkNode(cnode);
    //     cnode = temp;
    // }
    res.set(cnode);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Base of unary operations with signature: node * integer -> node
/////////////////////////////////////////////////////////////////////////////////////////

base_NItoN::base_NItoN(MEDDLY::opname* opcode, MEDDLY::expert_forest* _argF,
                       MEDDLY::expert_forest* _resF)
/**/ : MEDDLY::operation(opcode, 1), argF(_argF), resF(_resF)
{
    registerInForest(argF);
    registerInForest(resF);

    MEDDLY::ct_entry_type* et;
    et = new MEDDLY::ct_entry_type(opcode->getName(), "NI:N");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(3, resF);
    registerEntryType(0, et);
    buildCTs();

    mddUnion = MEDDLY::getOperation(MEDDLY::UNION, resF, resF, resF);  
}

base_NItoN::~base_NItoN() {
    unregisterInForest(argF);
    unregisterInForest(resF);   
}

////////////////////////////////////////////////////////////////////////////////////////

inline void 
base_NItoN::saveResult(MEDDLY::ct_entry_key* key, 
                       MEDDLY::node_handle a, const int b, MEDDLY::node_handle c) 
{
    CTresult[0].reset();
    CTresult[0].writeN(c);
    CT0->addEntry(key, CTresult[0]);
}

/////////////////////////////////////////////////////////////////////////////////////////

inline MEDDLY::ct_entry_key* 
base_NItoN::findResult(MEDDLY::node_handle a, const int b, MEDDLY::node_handle &c) 
{
    MEDDLY::ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    assert(CTsrch);
    CTsrch->writeN(a);
    CTsrch->writeI(b);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    c = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////

bool base_NItoN::checkForestCompatibility() const {
    return argF==resF;
}

/////////////////////////////////////////////////////////////////////////////////////////

void base_NItoN::computeDDEdge(const MEDDLY::dd_edge &a, const int b, MEDDLY::dd_edge &res)
{
    MEDDLY::node_handle cnode = compute(a.getNode(), b);
//   const int num_levels = resF->getDomain()->getNumVariables();
//   if ( userFlag && resF->isQuasiReduced() && cnode != resF->getTransparentNode()
//     && resF->getNodeLevel(cnode) < num_levels) {
//     node_handle temp = ((mt_forest*)resF)->makeNodeAtLevel(num_levels, cnode);
//     resF->unlinkNode(cnode);
//     cnode = temp;
//   }
    res.set(cnode);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Base class for unary operators performing: Node -> Node 
/////////////////////////////////////////////////////////////////////////////////////////

base_NtoN::base_NtoN(MEDDLY::unary_opname* opcode, MEDDLY::expert_forest* arg1, 
                     MEDDLY::expert_forest* res)
  : MEDDLY::unary_operation(opcode, 1, arg1, res)
{
    MEDDLY::ct_entry_type* et;
    et = new MEDDLY::ct_entry_type(opcode->getName(), "N:N");
    et->setForestForSlot(0, arg1);
    et->setForestForSlot(2, res);
    registerEntryType(0, et);
    buildCTs();

    mddUnion = MEDDLY::getOperation(MEDDLY::UNION, res, res, res);
}

/////////////////////////////////////////////////////////////////////////////////////////

inline MEDDLY::ct_entry_key* 
base_NtoN::findResult(MEDDLY::node_handle a, MEDDLY::node_handle &c) 
{
    MEDDLY::ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    assert(CTsrch);
    CTsrch->writeN(a);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    c = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
}

inline void 
base_NtoN::saveResult(MEDDLY::ct_entry_key* key, 
                      MEDDLY::node_handle a, MEDDLY::node_handle c) 
{
    CTresult[0].reset();
    CTresult[0].writeN(c);
    CT0->addEntry(key, CTresult[0]);
}

/////////////////////////////////////////////////////////////////////////////////////////

void base_NtoN::computeDDEdge(const MEDDLY::dd_edge &ar1, 
                              MEDDLY::dd_edge &res, bool userFlag)
{
    MEDDLY::node_handle cnode = compute(ar1.getNode());
    // const int num_levels = resF->getDomain()->getNumVariables();
    // if (userFlag && resF->isQuasiReduced() && cnode != resF->getTransparentNode()
    //     && resF->getNodeLevel(cnode) < num_levels) 
    // {
    //     MEDDLY::node_handle temp = (resF)->makeNodeAtLevel(num_levels, cnode);
    //     resF->unlinkNode(cnode);
    //     cnode = temp;
    // }
    res.set(cnode);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Base class for unary operators performing: Node -> int 
/////////////////////////////////////////////////////////////////////////////////////////

base_NtoI::base_NtoI(MEDDLY::unary_opname* opcode, MEDDLY::expert_forest* arg1)
  : MEDDLY::unary_operation(opcode, 1, arg1, MEDDLY::opnd_type::INTEGER)
{
    MEDDLY::ct_entry_type* et;
    et = new MEDDLY::ct_entry_type(opcode->getName(), "N:I");
    et->setForestForSlot(0, arg1);
    registerEntryType(0, et);
    buildCTs();
}

/////////////////////////////////////////////////////////////////////////////////////////

inline MEDDLY::ct_entry_key* 
base_NtoI::findResult(MEDDLY::node_handle a, int &c) 
{
    MEDDLY::ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    assert(CTsrch);
    CTsrch->writeN(a);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    c = CTresult[0].readI();
    CT0->recycle(CTsrch);
    return 0;
}

inline void 
base_NtoI::saveResult(MEDDLY::ct_entry_key* key, MEDDLY::node_handle a, int c) 
{
    CTresult[0].reset();
    CTresult[0].writeI(c);
    CT0->addEntry(key, CTresult[0]);
}

/////////////////////////////////////////////////////////////////////////////////////////

void base_NtoI::computeDDEdge(const MEDDLY::dd_edge &ar1, 
                              MEDDLY::dd_edge &res, bool userFlag)
{ throw; }

void base_NtoI::computeDDEdge(const MEDDLY::dd_edge &ar1, 
                              int &res, bool userFlag)
{
    res = compute(ar1.getNode());
}

/////////////////////////////////////////////////////////////////////////////////////////












/////////////////////////////////////////////////////////////////////////////////////////
// S-Vectors
/////////////////////////////////////////////////////////////////////////////////////////

s_vectors::s_vectors(MEDDLY::binary_opname* opcode, MEDDLY::expert_forest* arg1, 
                     MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res,
                     const s_vectors_table *tab, bool isConf, bool invertB, 
                     const sv_sign sign_of_sum, size_t lambda)
  : base_NNtoN(opcode, arg1, arg2, res), 
    p_table(tab), isPotentiallyConformant(isConf), invertB(invertB), 
    sign_of_sum(sign_of_sum), lambda(lambda)
{ 
    if (!invertB)
        operationCommutes();
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::node_handle s_vectors::compute(MEDDLY::node_handle a, MEDDLY::node_handle b) {
    if (a==0 || b==0) return 0;
    if (a==-1 && b==-1) return isPotentiallyConformant ? 0 : -1;

    MEDDLY::node_handle result;
    MEDDLY::ct_entry_key* key = findResult(a, b, result);
    if (nullptr==key)
        return result;

    const int a_level = arg1F->getNodeLevel(a);
    const int b_level = arg2F->getNodeLevel(b);
    assert(a_level == b_level);
    const int res_level = std::max(a_level, b_level);

    MEDDLY::unpacked_node *A = (a_level < res_level)
        ? MEDDLY::unpacked_node::newRedundant(arg1F, res_level, a, false)
        : arg1F->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    MEDDLY::unpacked_node *B = (b_level < res_level)
        ? MEDDLY::unpacked_node::newRedundant(arg2F, res_level, b, false)
        : arg2F->newUnpacked(b, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // if (a_level < res_level) 
    //     A->initRedundant(arg1F, res_level, a, false);
    // else
    //     arg1F->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);
    // MEDDLY::unpacked_node *B = MEDDLY::unpacked_node::New();
    // if (b_level < res_level) 
    //     B->initRedundant(arg1F, res_level, b, false);
    // else
    //     arg1F->unpackNode(B, b, MEDDLY::FULL_OR_SPARSE);
    // MEDDLY::unpacked_node *A = (a_level < res_level) 
    //     ? MEDDLY::unpacked_node::newRedundant(arg1F, res_level, a, false)
    //     : MEDDLY::unpacked_node::newFromNode(arg1F, a, false); 
    // MEDDLY::unpacked_node *B = (b_level < res_level)
    //     ? MEDDLY::unpacked_node::newRedundant(arg2F, res_level, b, false)
    //     : MEDDLY::unpacked_node::newFromNode(arg2F, b, false);

    const size_t a_size = get_node_size(A);
    const size_t b_size = get_node_size(B);
    const int res_size = a_size + b_size + 1; // with +1 because we are encoding integers
    check_level_bound(resF, res_level, res_size);
    
    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, res_level, res_size);

    const bool a_full = A->isFull(), b_full = B->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        for (size_t j = 0; j < (b_full ? b_size : B->getNNZs()); j++) { // for each b
            if (b_full && 0==B->d(j))
                continue;
            int b_val = NodeToZ(b_full ? j : B->i(j));

            if (invertB)
                b_val = -b_val;

            int ab_sum = add_exact(a_val, b_val); // a + b or a - b
            int ab_sign_prod = sign3(a_val) * sign3(b_val); // a * b
            // canonicalize the sum (if not yet decided)
            sv_sign curr_sign_of_sum = sign_of_sum;
            if (sign_of_sum == SVS_UNDECIDED && 0!=ab_sum) { 
                curr_sign_of_sum = (ab_sum > 0) ? SVS_POS : SVS_NEG;
            }
            // compute -(a+b) or -(a-b)
            if (curr_sign_of_sum == SVS_NEG) {
                ab_sum = -ab_sum;
                // ab_prod = -ab_prod;
            }

            int ab_sum_idx = ZtoNode(ab_sum);
            if (ab_sum_idx >= res_size)
                cout << "ab_sum_idx:"<<ab_sum_idx<<" res_size:"<<res_size<<endl;
            assert(ab_sum_idx < res_size);

            bool ij_conf = (ab_sign_prod >= 0);
            // if (curr_sign_of_sum == SVS_NEG)
            //     ij_conf = -(multiply_exact(a_val, b_val)) >= 0;
            // else
            //     ij_conf = (multiply_exact(a_val, b_val)) >= 0;



            bool down_is_conf = isPotentiallyConformant && ij_conf;

            bool do_sum = true;
            if (lambda==0) {
            }
            else if (p_table->pivot_order->is_same_as_lambda(lambda, res_level)) {
                do_sum = !ij_conf;
            }
            else if (p_table->pivot_order->is_below_lambda(lambda, res_level)) {
                do_sum = ij_conf;
            }
            else if (p_table->pivot_order->is_above_lambda(lambda, res_level)) {
            }
            // cout << "a_level:"<<a_level<<" a_val:"<<a_val<<" b_val:"<<b_val<<" ij_conf:"<<ij_conf<<" do_sum:"<<do_sum<<endl;

            if (do_sum) {
                MEDDLY::node_handle sum;
                s_vectors *next_op = p_table->get_op(lambda, down_is_conf, invertB, curr_sign_of_sum);
                sum = next_op->compute(A->d(i), B->d(j));
                unionNodes(C, sum, ab_sum_idx, resF, mddUnion);
            }
        }
    }

    // cleanup
    MEDDLY::unpacked_node::recycle(B);
    MEDDLY::unpacked_node::recycle(A);

    // reduce and return result
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, b, result);
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

s_vectors* 
s_vectors_opname::buildOperation(MEDDLY::expert_forest* a1, 
        MEDDLY::expert_forest* a2, MEDDLY::expert_forest* r,
        const s_vectors_table *tab, bool isConf, bool invertB, 
        const sv_sign sign_of_sum, size_t lambda)
{
    if (0==a1 || 0==a2 || 0==r) return 0;

    if ((a1->getDomain() != r->getDomain()) || (a2->getDomain() != r->getDomain()))
        throw MEDDLY::error(MEDDLY::error::DOMAIN_MISMATCH, __FILE__, __LINE__);

    if (a1->isForRelations() || r->isForRelations() || a2->isForRelations() ||
        (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
        (a2->getEdgeLabeling() != r->getEdgeLabeling()))
        throw MEDDLY::error(MEDDLY::error::TYPE_MISMATCH, __FILE__, __LINE__);

    if (r->getEdgeLabeling() == MEDDLY::edge_labeling::MULTI_TERMINAL) {
        if (r->isForRelations())
            throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED);
        return new s_vectors(this, a1, a2, r, tab, isConf, invertB, sign_of_sum, lambda);
    }
    throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

MEDDLY::binary_operation* 
s_vectors_opname::buildOperation(MEDDLY::expert_forest* arg1,
                                 MEDDLY::expert_forest* arg2, 
                                 MEDDLY::expert_forest* res) 
{
    throw std::exception(); // Unimplemented
}

/////////////////////////////////////////////////////////////////////////////////////////

s_vectors_table::s_vectors_table(MEDDLY::expert_forest* forest, 
                                 const variable_order *pivot_order) 
/**/ : forest(forest), pivot_order(pivot_order)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

s_vectors* 
s_vectors_table::get_op(size_t level, bool isPotentiallyConformant, 
                        bool invertB, const sv_sign sign_of_sum) const
{
    if (table.empty())
        table.resize(forest->getNumVariables() + 1);
    if (table[level].empty())
        table[level].resize(2);
    if (table[level][isPotentiallyConformant].empty())
        table[level][isPotentiallyConformant].resize(2);
    if (table[level][isPotentiallyConformant][invertB].empty())
        table[level][isPotentiallyConformant][invertB].resize(SVS_TOTAL);
    if (table[level][isPotentiallyConformant][invertB][sign_of_sum] == nullptr)
        table[level][isPotentiallyConformant][invertB][sign_of_sum] =
            S_VECTORS_OPNAME->buildOperation(forest, forest, forest, this, 
                                             isPotentiallyConformant, invertB,
                                             sign_of_sum, level);
    
    return table[level][isPotentiallyConformant][invertB][sign_of_sum];
}

/////////////////////////////////////////////////////////////////////////////////////////









/////////////////////////////////////////////////////////////////////////////////////////
// Less Equal Squared Operator
// Get the elements using the less-equal-but-not-equal-squared operator
/////////////////////////////////////////////////////////////////////////////////////////

lesseq_sq::lesseq_sq(MEDDLY::binary_opname* opcode, 
/**/        MEDDLY::expert_forest* arg1, MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res,
/**/        const lesseq_sq_table* tab,
/**/        bool isEq, bool isB0,
/**/        bool _subtract, size_t lambda)
  : base_NNtoN(opcode, arg1, arg2, res), p_table(tab),
    isPotentiallyEqual(isEq), isBZero(isB0),
    subtract(_subtract), lambda(lambda)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::node_handle lesseq_sq::compute(MEDDLY::node_handle a, MEDDLY::node_handle b) {
    if (a==0) return 0;
    if (b==0) return 0;
    if (a==-1) return isPotentiallyEqual || isBZero ? 0 : -1;
    assert(b != -1);

    MEDDLY::node_handle result;
    MEDDLY::ct_entry_key* key = findResult(a, b, result);
    if (nullptr==key)
        return result;

    const int a_level = arg1F->getNodeLevel(a);
    const int b_level = arg2F->getNodeLevel(b);
    assert(a_level == b_level);
    const int res_level = std::max(a_level, b_level);

    MEDDLY::unpacked_node *A = (a_level < res_level)
        ? MEDDLY::unpacked_node::newRedundant(arg1F, res_level, a, false)
        : arg1F->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    MEDDLY::unpacked_node *B = (b_level < res_level)
        ? MEDDLY::unpacked_node::newRedundant(arg2F, res_level, b, false)
        : arg2F->newUnpacked(b, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // if (a_level < res_level) 
    //     A->initRedundant(arg1F, res_level, a, false);
    // else
    //     arg1F->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);
    // MEDDLY::unpacked_node *B = MEDDLY::unpacked_node::New();
    // if (b_level < res_level) 
    //     B->initRedundant(arg1F, res_level, b, false);
    // else
    //     arg1F->unpackNode(B, b, MEDDLY::FULL_OR_SPARSE);
    // MEDDLY::unpacked_node *A = (a_level < res_level) 
    //     ? MEDDLY::unpacked_node::newRedundant(arg1F, res_level, a, false)
    //     : MEDDLY::unpacked_node::newFromNode(arg1F, a, false);
    // MEDDLY::unpacked_node *B = (b_level < res_level)
    //     ? MEDDLY::unpacked_node::newRedundant(arg2F, res_level, b, false)
    //     : MEDDLY::unpacked_node::newFromNode(arg2F, b, false);

    const size_t a_size = get_node_size(A);
    const size_t b_size = get_node_size(B);
    // const size_t res_size = a_size;
    size_t res_size;
    if (lambda != 0 && 
        p_table->pivot_order->is_above_lambda(lambda, res_level) && subtract) 
    {
        res_size = a_size + b_size;
    }
    else res_size = a_size;
    check_level_bound(resF, res_level, res_size);

    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, res_level, res_size);

    const bool a_full = A->isFull(), b_full = B->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        for (size_t j = 0; j < (b_full ? b_size : B->getNNZs()); j++) { // for each b
            if (b_full && 0==B->d(j))
                continue;
            int b_val = NodeToZ(b_full ? j : B->i(j));

            bool is_lesseq_sq;
            bool is_pot_eq = isPotentiallyEqual;
            bool is_b_pot_zero = isBZero;
            if (lambda != 0 && p_table->pivot_order->is_above_lambda(lambda, res_level)) {
                is_lesseq_sq = true;
            }
            else {
            //     // check that i <= j and both are conformal
                int ab_sign_prod = sign3(a_val) * sign3(b_val);
                is_lesseq_sq = (abs(b_val) <= abs(a_val) && ab_sign_prod >= 0);
                is_pot_eq &= (a_val == b_val);
                is_b_pot_zero &= (0 == b_val);
            }

            if (is_lesseq_sq) 
            { 
                int idx = ZtoNode(subtract ? subtract_exact(a_val, b_val) : a_val);

                // Determine next operator in chain
                bool down_is_eq = isPotentiallyEqual && is_pot_eq;
                bool down_is_b0 = isBZero && is_b_pot_zero;
                lesseq_sq* next_op = p_table->get_op(lambda, down_is_eq, 
                                                     down_is_b0, subtract);

                MEDDLY::node_handle leq_ij = next_op->compute(A->d(i), B->d(j));
                unionNodes(C, leq_ij, idx, resF, mddUnion);
            }
        }
    }

    // cleanup
    MEDDLY::unpacked_node::recycle(B);
    MEDDLY::unpacked_node::recycle(A);
    // reduce and return result
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, b, result);
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

lesseq_sq* 
lesseq_sq_opname::buildOperation(MEDDLY::expert_forest* a1, 
                                          MEDDLY::expert_forest* a2, 
                                          MEDDLY::expert_forest* r,
                                          const lesseq_sq_table* tab,
                                          bool isEq, bool isB0, bool subtract,
                                          size_t lambda)
{
    if (0==a1 || 0==a2 || 0==r) return 0;

    if ((a1->getDomain() != r->getDomain()) || (a2->getDomain() != r->getDomain()))
        throw MEDDLY::error(MEDDLY::error::DOMAIN_MISMATCH, __FILE__, __LINE__);

    if (a1->isForRelations() || r->isForRelations() || a2->isForRelations() ||
        (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
        (a2->getEdgeLabeling() != r->getEdgeLabeling()))
        throw MEDDLY::error(MEDDLY::error::TYPE_MISMATCH, __FILE__, __LINE__);

    if (r->getEdgeLabeling() == MEDDLY::edge_labeling::MULTI_TERMINAL) {
        if (r->isForRelations())
            throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED);
        return new lesseq_sq(this, a1, a2, r, tab, isEq, isB0,
                             subtract, lambda);
    }
    throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

MEDDLY::binary_operation* 
lesseq_sq_opname::buildOperation(MEDDLY::expert_forest* a1, 
                                          MEDDLY::expert_forest* a2, 
                                          MEDDLY::expert_forest* r)
{
    throw std::exception(); // Unimplemented
}

/////////////////////////////////////////////////////////////////////////////////////////

lesseq_sq_table::lesseq_sq_table(MEDDLY::expert_forest* forest,
                                 const variable_order *pivot_order) 
/**/ : forest(forest), pivot_order(pivot_order)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

lesseq_sq* 
lesseq_sq_table:: get_op(size_t level, bool isPotentiallyEqual, 
                         bool isBZero, bool subtract) const 
{
    if (table.empty())
        table.resize(forest->getNumVariables() + 1);
    if (table[level].empty())
        table[level].resize(2);
    if (table[level][isPotentiallyEqual].empty())
        table[level][isPotentiallyEqual].resize(2);
    if (table[level][isPotentiallyEqual][isBZero].empty())
        table[level][isPotentiallyEqual][isBZero].resize(2);
    if (table[level][isPotentiallyEqual][isBZero][subtract] == nullptr)
        table[level][isPotentiallyEqual][isBZero][subtract] = 
            LESSEQ_SQ_OPNAME->buildOperation(forest, forest, forest, this, 
                                             isPotentiallyEqual, isBZero, 
                                             subtract, level);

    return table[level][isPotentiallyEqual][isBZero][subtract];
}

/////////////////////////////////////////////////////////////////////////////////////////


















// /////////////////////////////////////////////////////////////////////////////////////////
// // Completion Procedure
// // Remove elements that are equal between [0..j-1] and smaller at var j
// /////////////////////////////////////////////////////////////////////////////////////////

// compl_proc::compl_proc(MEDDLY::binary_opname* opcode, 
// /**/                   MEDDLY::expert_forest* arg1, MEDDLY::expert_forest* arg2, 
// /**/                   MEDDLY::expert_forest* res,
// /**/                   const compl_proc_table* tab,
// /**/                   size_t _rlvl)
//   : base_NNtoN(opcode, arg1, arg2, res), p_table(tab), restricted_level(_rlvl)
// { }

// /////////////////////////////////////////////////////////////////////////////////////////

// MEDDLY::node_handle compl_proc::compute(MEDDLY::node_handle a, MEDDLY::node_handle b) {
//     if (a==0 || b==0) return 0;
//     if (a==-1 && b==-1) return -1;

//     MEDDLY::node_handle result;
//     MEDDLY::ct_entry_key* key = findResult(a, b, result);
//     if (nullptr==key)
//         return result;

//     const int a_level = arg1F->getNodeLevel(a);
//     const int b_level = arg2F->getNodeLevel(b);
//     assert(a_level == b_level);
//     const int res_level = std::max(a_level, b_level);

//     MEDDLY::unpacked_node *A = (a_level < res_level) 
//         ? MEDDLY::unpacked_node::newRedundant(arg1F, res_level, a, false)
//         : MEDDLY::unpacked_node::newFromNode(arg1F, a, false);
//     MEDDLY::unpacked_node *B = (b_level < res_level)
//         ? MEDDLY::unpacked_node::newRedundant(arg2F, res_level, b, false)
//         : MEDDLY::unpacked_node::newFromNode(arg2F, b, false);

//     const size_t a_size = get_node_size(A);
//     const size_t b_size = get_node_size(B);
//     // const size_t res_size = a_size;
//     size_t res_size = a_size + b_size;
//     check_level_bound(resF, res_level, res_size);

//     MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, res_level, res_size);

//     const bool a_full = A->isFull(), b_full = B->isFull();

//     for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
//         if (a_full && 0==A->d(i))
//             continue;
//         int a_val = NodeToZ(a_full ? i : A->i(i));

//         for (size_t j = 0; j < (b_full ? b_size : B->getNNZs()); j++) { // for each b
//             if (b_full && 0==B->d(j))
//                 continue;
//             int b_val = NodeToZ(b_full ? j : B->i(j));

//             bool is_selected_for_removal;
//             if (res_level > restricted_level) { // take everything
//                 is_selected_for_removal = true;
//             }
//             else if (res_level == restricted_level) { // check if greater in absolute value
//                 is_selected_for_removal = (abs(a_val) > abs(b_val)) && (a_val*b_val > 0);
//             }
//             else {  // check equality
//                 is_selected_for_removal = (a_val == b_val);
//             }
//             // bool is_compl_proc;
//             // bool is_pot_eq = isPotentiallyEqual;
//             // bool is_b_pot_zero = isBZero;
//             // if (restricted_level != 0 && res_level > restricted_level) {
//             //     is_compl_proc = true;
//             // }
//             // else {
//             // //     // check that i <= j and both are conformal
//             //     is_compl_proc = (abs(b_val) <= abs(a_val) && multiply_exact(a_val, b_val) >= 0);
//             //     is_pot_eq &= (a_val == b_val);
//             //     is_b_pot_zero &= (0 == b_val);
//             // }

//             if (is_selected_for_removal) 
//             { 
//                 int idx = ZtoNode(a_val);

//                 // // Determine next operator in chain
//                 // bool down_is_eq = isPotentiallyEqual && is_pot_eq;
//                 // bool down_is_b0 = isBZero && is_b_pot_zero;
//                 // // bool down_level = (restricted_level!=0 && res_level > restricted_level) ? restricted_level : 0;
//                 // // bool down_to_neq = isPotentiallyEqual && !is_pot_eq;
//                 // // bool down_to_nz = isBZero && !is_b_pot_zero;
//                 // compl_proc* next_op = p_table->get_op(restricted_level, down_is_eq, 
//                 //                                      down_is_b0, subtract);

//                 MEDDLY::node_handle res_ij = compute(A->d(i), B->d(j));
//                 unionNodes(C, res_ij, idx, resF, mddUnion);
//             }
//         }
//     }

//     // cleanup
//     MEDDLY::unpacked_node::recycle(B);
//     MEDDLY::unpacked_node::recycle(A);
//     // reduce and return result
//     result = resF->createReducedNode(-1, C);
//     saveResult(key, a, b, result);
//     return result;
// }

// /////////////////////////////////////////////////////////////////////////////////////////

// compl_proc* 
// compl_proc_opname::buildOperation(MEDDLY::expert_forest* a1, 
//                                   MEDDLY::expert_forest* a2, 
//                                   MEDDLY::expert_forest* r,
//                                   const compl_proc_table* tab,
//                                   size_t restricted_level)
// {
//     if (0==a1 || 0==a2 || 0==r) return 0;

//     if ((a1->getDomain() != r->getDomain()) || (a2->getDomain() != r->getDomain()))
//         throw MEDDLY::error(MEDDLY::error::DOMAIN_MISMATCH, __FILE__, __LINE__);

//     if (a1->isForRelations() || r->isForRelations() || a2->isForRelations() ||
//         (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
//         (a2->getEdgeLabeling() != r->getEdgeLabeling()))
//         throw MEDDLY::error(MEDDLY::error::TYPE_MISMATCH, __FILE__, __LINE__);

//     if (r->getEdgeLabeling() == MEDDLY::edge_labeling::MULTI_TERMINAL) {
//         if (r->isForRelations())
//             throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED);
//         return new compl_proc(this, a1, a2, r, tab, restricted_level);
//     }
//     throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED, __FILE__, __LINE__);
// }

// MEDDLY::binary_operation* 
// compl_proc_opname::buildOperation(MEDDLY::expert_forest* a1, 
//                                           MEDDLY::expert_forest* a2, 
//                                           MEDDLY::expert_forest* r)
// {
//     throw std::exception(); // Unimplemented
// }

// /////////////////////////////////////////////////////////////////////////////////////////

// compl_proc_table::compl_proc_table(MEDDLY::expert_forest* forest) {
//     const size_t num_vars = forest->getNumVariables();
//     table.resize(num_vars + 1);
//     for (size_t lvl=0; lvl<=num_vars; ++lvl) {
//         table[lvl] =
//             COMPL_PROC_OPNAME->buildOperation(forest, forest, forest,
//                                               this, lvl);
//     }
// }

// /////////////////////////////////////////////////////////////////////////////////////////






















/////////////////////////////////////////////////////////////////////////////////////////
// Canonicalize the sign of a MDD
/////////////////////////////////////////////////////////////////////////////////////////

sign_canon_mdd_opname::sign_canon_mdd_opname() 
/**/ : unary_opname("CanonicalizeSign")
{ }

sign_canon_mdd_op* 
sign_canon_mdd_opname::buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF, 
                                      const sign_canon_table* p_table, bool isZero)
{
    if (0==arF || 0==resF) return 0;
    
    return new sign_canon_mdd_op(this, (MEDDLY::expert_forest*)arF, 
                                 (MEDDLY::expert_forest*)resF, p_table, isZero);
}

sign_canon_mdd_op::~sign_canon_mdd_op() { }

/////////////////////////////////////////////////////////////////////////////////////////

sign_canon_mdd_op::sign_canon_mdd_op(MEDDLY::unary_opname* opcode, MEDDLY::expert_forest* _argF,
                                     MEDDLY::expert_forest* _resF, const sign_canon_table* p_table, bool isZero)
/**/ : base_NtoN(opcode, _argF, _resF), p_table(p_table), isZero(isZero)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::node_handle sign_canon_mdd_op::compute(MEDDLY::node_handle a) 
{
    if (a==0 || a==-1) return a;

    MEDDLY::node_handle result;
    MEDDLY::ct_entry_key* CTsrch;
    MEDDLY::ct_entry_key* key = findResult(a, result);
    if (nullptr==key)
        return result;

    // Read the node and swap the sign of the negative values        
    MEDDLY::unpacked_node* A = argF->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node* A = MEDDLY::unpacked_node::newFromNode(argF, a, false);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // argF->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);
    const int a_level = argF->getNodeLevel(a);
    const size_t a_size = get_node_size(A);

    const int res_size = max(a_size, (size_t)ZtoNode(-NodeToZ(a_size-1)) + 1);
    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, a_level, res_size);

    const bool a_full = A->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        if (isZero == true) { // following a 0-path, not found yet the first nnz
            if (a_val > 0) { // keep as-is, end recursive descent
                unionNodes(C, resF->linkNode(A->d(i)), 
                           ZtoNode(a_val), resF, mddUnion);
            }
            else if (a_val == 0) { // continue descending
                unionNodes(C, compute(A->d(i)), 
                           ZtoNode(a_val), resF, mddUnion);
            }
            else { // invert vectors
                unionNodes(C, p_table->get_op(false)->compute(A->d(i)), 
                           ZtoNode(-a_val), resF, mddUnion);
            }
        }
        else { // keep inverting the vector
            unionNodes(C, compute(A->d(i)), 
                       ZtoNode(-a_val), resF, mddUnion);
        }
    }

    // // Cleanup
    MEDDLY::unpacked_node::recycle(A);
    // Save in cache
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, result);
    return result;
 }

/////////////////////////////////////////////////////////////////////////////////////////

sign_canon_table::sign_canon_table(MEDDLY::expert_forest* forest) 
/**/ : forest(forest)
{ }

sign_canon_mdd_op* sign_canon_table::get_op(bool isZero) const
{
    if (table.empty())
        table.resize(2);
    if (table[isZero] == nullptr)
        table[isZero] = SIGN_CANON_OPNAME->buildOperation(forest, forest, this, isZero);

    return table[isZero];
}

/////////////////////////////////////////////////////////////////////////////////////////














/////////////////////////////////////////////////////////////////////////////////////////
// Multiply vectors by a constant
/////////////////////////////////////////////////////////////////////////////////////////

vmult_opname::vmult_opname() 
/**/ : unary_opname("VectorMultiply")
{ }

MEDDLY::operation* 
vmult_opname::buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF)
{
    if (0==arF || 0==resF) return 0;
    
    return new vmult_op(this, (MEDDLY::expert_forest*)arF, 
                             (MEDDLY::expert_forest*)resF);
}

vmult_op::~vmult_op() { }

/////////////////////////////////////////////////////////////////////////////////////////

vmult_op::vmult_op(MEDDLY::opname* opcode, MEDDLY::expert_forest* _argF,
                   MEDDLY::expert_forest* _resF)
/**/ : base_NItoN(opcode, _argF, _resF)
{ }


/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::node_handle vmult_op::compute(MEDDLY::node_handle a, const int multiplier) 
{
    if (a==0 || a==-1) return a;
    if (multiplier==1) return resF->linkNode(a);

    MEDDLY::node_handle result;
    MEDDLY::ct_entry_key* CTsrch;
    MEDDLY::ct_entry_key* key = findResult(a, multiplier, result);
    if (nullptr==key)
        return result;

    // // Read the node and accumulate the LCM of the GCDs
    // MEDDLY::unpacked_node* A = MEDDLY::unpacked_node::newFromNode(argF, a, false);
    MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    argF->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);
    const int a_level = argF->getNodeLevel(a);
    const size_t a_size = get_node_size(A);

    const int res_size = a_size * abs(multiplier) + (multiplier<0 ? 1 : 0);
    check_level_bound(resF, a_level, res_size);
    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, a_level, res_size);

    const bool a_full = A->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        int r_val = ZtoNode(multiply_exact(multiplier, a_val));

        unionNodes(C, compute(A->d(i), multiplier), 
                   r_val, resF, mddUnion);
    }

    // // Cleanup
    MEDDLY::unpacked_node::recycle(A);
    // Save in cache
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, multiplier, result);
    return result;
 }

/////////////////////////////////////////////////////////////////////////////////////////




















/////////////////////////////////////////////////////////////////////////////////////////
// Canonicalize vectors
// Divide all paths that are divisible by a given number
/////////////////////////////////////////////////////////////////////////////////////////

vcanon_mdd_opname::vcanon_mdd_opname(bool _divide) 
/**/ : unary_opname("CanonicalizeGCD"), divide(_divide)
{ }

MEDDLY::operation* 
vcanon_mdd_opname::buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF)
{
    if (0==arF || 0==resF) return 0;
    
    return new vcanon_mdd_op(this, (MEDDLY::expert_forest*)arF, 
                             (MEDDLY::expert_forest*)resF, divide);
}

vcanon_mdd_op::~vcanon_mdd_op() { }

/////////////////////////////////////////////////////////////////////////////////////////

vcanon_mdd_op::vcanon_mdd_op(MEDDLY::opname* opcode, MEDDLY::expert_forest* _argF,
                             MEDDLY::expert_forest* _resF, bool _divide)
/**/ : base_NItoN(opcode, _argF, _resF), divide(_divide)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::node_handle vcanon_mdd_op::compute(MEDDLY::node_handle a, const int divisor) 
{
    if (a==0 || a==-1) return a;
    if (divisor==1) return resF->linkNode(a);

    MEDDLY::node_handle result;
    MEDDLY::ct_entry_key* CTsrch;
    MEDDLY::ct_entry_key* key = findResult(a, divisor, result);
    if (nullptr==key)
        return result;

    // // Read the node and accumulate the LCM of the GCDs
    MEDDLY::unpacked_node* A = argF->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node* A = MEDDLY::unpacked_node::newFromNode(argF, a, false);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // argF->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);
    const int a_level = argF->getNodeLevel(a);
    const size_t a_size = get_node_size(A);

    const int res_size = a_size;
    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, a_level, res_size);

    const bool a_full = A->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        if (0 == (abs(a_val) % divisor)) // path is divisible by @divisor
            unionNodes(C, compute(A->d(i), divisor), 
                        divide ? ZtoNode(a_val / divisor) : ZtoNode(a_val), 
                        resF, mddUnion);

    }

    // // Cleanup
    MEDDLY::unpacked_node::recycle(A);
    // Save in cache
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, divisor, result);
    return result;
 }

/////////////////////////////////////////////////////////////////////////////////////////





























/////////////////////////////////////////////////////////////////////////////////////////
// Divisors of DD paths
/////////////////////////////////////////////////////////////////////////////////////////

divisors_finder_mdd_opname::divisors_finder_mdd_opname() 
/**/ : unary_opname("DivisorsFinderMDD") { }

MEDDLY::unary_operation* 
divisors_finder_mdd_opname::buildOperation(MEDDLY::forest* arg)
{
    if (0==arg) return 0;
    
    return new divisors_finder_mdd_op(this, (MEDDLY::expert_forest*)arg);
}

/////////////////////////////////////////////////////////////////////////////////////////

divisors_finder_mdd_op::divisors_finder_mdd_op(MEDDLY::unary_opname* oc, MEDDLY::expert_forest* arg) 
/**/ : base_NtoI(oc, arg)
{
    lut.insert(std::vector<int>{});  // 0 terminal
    lut.insert(std::vector<int>{0}); // 1 terminal
}

/////////////////////////////////////////////////////////////////////////////////////////

divisors_finder_mdd_op::lut_key 
divisors_finder_mdd_op::compute(MEDDLY::node_handle a)
{
    if (0 == a) 
        return 0; // return { }
    if (-1 == a) 
        return 1; // return { 0 }
    
    lut_key res;
    MEDDLY::ct_entry_key* CTsrch;
    MEDDLY::ct_entry_key* key = findResult(a, res);
    if (nullptr==key)
        return res;

    // Read the node and accumulate the set of all GCDs
    std::vector<int> divisors;
    MEDDLY::unpacked_node* A = argF->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node* A = MEDDLY::unpacked_node::newFromNode(argF, a, false);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // argF->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);

    const size_t a_size = get_node_size(A);
    const bool a_full = A->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        lut_key lk_i = compute(A->d(i));
        const std::vector<int>& div_i = lut.look_up(lk_i);
        for (int v : div_i) {
            int new_div_i = gcd(v, abs(a_val));
            // cout << "gcd("<<v<<", "<<abs(a_val)<<")="<<new_div_i<<" ";
            if ((new_div_i > 1) || (a_val==0 && new_div_i==0))
                divisors.push_back(new_div_i);
        }
        // cout << endl;
    }

    // Sort and make the list of divisors unique
    std::sort(divisors.begin(), divisors.end(), [](int i, int j){ return i>j;});
    divisors.erase(std::unique(divisors.begin(), divisors.end()), divisors.end());

    // Cleanup
    MEDDLY::unpacked_node::recycle(A);
    // Save in cache
    res = lut.insert(std::move(divisors));
    saveResult(key, a, res);
    return res;
}

/////////////////////////////////////////////////////////////////////////////////////////























/////////////////////////////////////////////////////////////////////////////////////////
// Less Equal Squared Operator
// Get the elements using the less-equal-but-not-equal-squared operator
/////////////////////////////////////////////////////////////////////////////////////////

support_inclusion_op::support_inclusion_op(MEDDLY::binary_opname* opcode, 
/**/        MEDDLY::expert_forest* arg1, MEDDLY::expert_forest* arg2, MEDDLY::expert_forest* res,
/**/        const support_inclusion_table* p_table, bool is_pot_eq_supp, bool subtract, size_t lambda)
  : base_NNtoN(opcode, arg1, arg2, res), p_table(p_table), 
    is_pot_eq_supp(is_pot_eq_supp), subtract(subtract), lambda(lambda)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::node_handle support_inclusion_op::compute(MEDDLY::node_handle a, MEDDLY::node_handle b) {
    if (a==0) return 0;
    if (b==0) return 0;
    if (a==-1) return is_pot_eq_supp ? 0 : -1;
    assert(b != -1);

    MEDDLY::node_handle result;
    MEDDLY::ct_entry_key* key = findResult(a, b, result);
    if (nullptr==key)
        return result;

    const int a_level = arg1F->getNodeLevel(a);
    const int b_level = arg2F->getNodeLevel(b);
    assert(a_level == b_level);
    const int res_level = std::max(a_level, b_level);

    MEDDLY::unpacked_node *A = (a_level < res_level)
        ? MEDDLY::unpacked_node::newRedundant(arg1F, res_level, a, false)
        : arg1F->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    MEDDLY::unpacked_node *B = (b_level < res_level)
        ? MEDDLY::unpacked_node::newRedundant(arg2F, res_level, b, false)
        : arg2F->newUnpacked(b, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // if (a_level < res_level) 
    //     A->initRedundant(arg1F, res_level, a, false);
    // else
    //     arg1F->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);
    // MEDDLY::unpacked_node *B = MEDDLY::unpacked_node::New();
    // if (b_level < res_level) 
    //     B->initRedundant(arg1F, res_level, b, false);
    // else
    //     arg1F->unpackNode(B, b, MEDDLY::FULL_OR_SPARSE);
    // MEDDLY::unpacked_node *A = (a_level < res_level) 
    //     ? MEDDLY::unpacked_node::newRedundant(arg1F, res_level, a, false)
    //     : MEDDLY::unpacked_node::newFromNode(arg1F, a, false);
    // MEDDLY::unpacked_node *B = (b_level < res_level)
    //     ? MEDDLY::unpacked_node::newRedundant(arg2F, res_level, b, false)
    //     : MEDDLY::unpacked_node::newFromNode(arg2F, b, false);

    const size_t a_size = get_node_size(A);
    const size_t b_size = get_node_size(B);
    const int res_size = subtract ? a_size + b_size : a_size;
    check_level_bound(resF, res_level, res_size);

    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, res_level, res_size);

    const bool a_full = A->isFull(), b_full = B->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        for (size_t j = 0; j < (b_full ? b_size : B->getNNZs()); j++) { // for each b
            if (b_full && 0==B->d(j))
                continue;
            int b_val = NodeToZ(b_full ? j : B->i(j));

            MEDDLY::node_handle n = 0; // n = a if supp(a) subseteq_neq supp(b)
            if (lambda != 0 && p_table->pivot_order->is_above_lambda(lambda, res_level)) {
                n = compute(A->d(i), B->d(j)); // no check
            }
            else if (a_val == 0) { // a==0 -> b can only be 0
                if (b_val == 0) {
                    n = compute(A->d(i), B->d(j));
                }
            }
            else { // a_val != 0 -> b can be any value or 0
                if (is_pot_eq_supp && b_val==0) // smaller support
                    n = p_table->get_op(lambda, false, subtract)->compute(A->d(i), B->d(j));
                else //if (multiply_exact(a_val, b_val) >= 0)
                    n = compute(A->d(i), B->d(j));
            }
            if (n != 0)
                unionNodes(C, n, ZtoNode(subtract ? a_val - b_val : a_val), 
                           resF, mddUnion);
        }
    }

    // cleanup
    MEDDLY::unpacked_node::recycle(B);
    MEDDLY::unpacked_node::recycle(A);
    // reduce and return result
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, b, result);
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

support_inclusion_op* 
support_inclusion_opname::buildOperation(MEDDLY::expert_forest* a1, 
                                         MEDDLY::expert_forest* a2, 
                                         MEDDLY::expert_forest* r,
                                         const support_inclusion_table* p_table,
                                         bool is_pot_eq_supp, bool subtract,
                                         size_t lambda)
{
    if (0==a1 || 0==a2 || 0==r) return 0;

    if ((a1->getDomain() != r->getDomain()) || (a2->getDomain() != r->getDomain()))
        throw MEDDLY::error(MEDDLY::error::DOMAIN_MISMATCH, __FILE__, __LINE__);

    if (a1->isForRelations() || r->isForRelations() || a2->isForRelations() ||
        (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
        (a2->getEdgeLabeling() != r->getEdgeLabeling()))
        throw MEDDLY::error(MEDDLY::error::TYPE_MISMATCH, __FILE__, __LINE__);

    if (r->getEdgeLabeling() == MEDDLY::edge_labeling::MULTI_TERMINAL) {
        if (r->isForRelations())
            throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED);
        return new support_inclusion_op(this, a1, a2, r, p_table, is_pot_eq_supp, subtract, lambda);
    }
    throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

MEDDLY::binary_operation* 
support_inclusion_opname::buildOperation(MEDDLY::expert_forest* a1, 
        MEDDLY::expert_forest* a2, MEDDLY::expert_forest* r)
{
    throw std::exception(); // Unimplemented
}

/////////////////////////////////////////////////////////////////////////////////////////

support_inclusion_table::support_inclusion_table(MEDDLY::expert_forest* forest,
                                 const variable_order *pivot_order) 
/**/ : forest(forest), pivot_order(pivot_order)
{ }

support_inclusion_op* 
support_inclusion_table::get_op(size_t level, bool isPotentiallyEqual, bool subtract) const
{
    if (table.empty())
        table.resize(forest->getNumVariables() + 1);
    if (table[level].empty())
        table[level].resize(2);
    if (table[level][isPotentiallyEqual].empty())
        table[level][isPotentiallyEqual].resize(2);
    if (table[level][isPotentiallyEqual][subtract] == nullptr)
        table[level][isPotentiallyEqual][subtract] = 
            SUPPORT_INCL_OPNAME->buildOperation(forest, forest, forest, this, 
                                                isPotentiallyEqual, subtract, level);

    return table[level][isPotentiallyEqual][subtract];

}

/////////////////////////////////////////////////////////////////////////////////////////





















/////////////////////////////////////////////////////////////////////////////////////////
// Level/Value selector
// Selector for a subset of an MDD having a given value at a specified level lambda
/////////////////////////////////////////////////////////////////////////////////////////

lv_selector_opname::lv_selector_opname() 
/**/ : unary_opname("LevelValueSelector")
{ }

lv_selector_op* 
lv_selector_opname::buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF, const size_t lambda)
{
    if (0==arF || 0==resF) return 0;
    
    return new lv_selector_op(this, (MEDDLY::expert_forest*)arF, 
                              (MEDDLY::expert_forest*)resF, lambda);
}

lv_selector_op::~lv_selector_op() { }

/////////////////////////////////////////////////////////////////////////////////////////

lv_selector_op::lv_selector_op(MEDDLY::opname* opcode, MEDDLY::expert_forest* _argF,
                             MEDDLY::expert_forest* _resF, const size_t lambda)
/**/ : base_NItoN(opcode, _argF, _resF), lambda(lambda)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::node_handle lv_selector_op::compute(MEDDLY::node_handle a, const int value) 
{
    if (a==0 || a==-1) return a;

    const int a_level = argF->getNodeLevel(a);
    if (a_level < lambda)
        return resF->linkNode(a);

    MEDDLY::node_handle result;
    MEDDLY::ct_entry_key* CTsrch;
    MEDDLY::ct_entry_key* key = findResult(a, value, result);
    if (nullptr==key)
        return result;

    // // Read the node and accumulate the LCM of the GCDs
    MEDDLY::unpacked_node* A = argF->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node* A = MEDDLY::unpacked_node::newFromNode(argF, a, false);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // argF->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);
    const size_t a_size = get_node_size(A);

    const int res_size = a_size;
    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, a_level, res_size);

    const bool a_full = A->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        // cout << "a_val:"<<a_val<<" a_level:"<<a_level<<" value:"<<value<<" lambda:"<<lambda<<endl;
        if (a_level != lambda || a_val == value) // select this path that has @value at level @lambda
            unionNodes(C, compute(A->d(i), value), ZtoNode(a_val), resF, mddUnion);
    }

    // Cleanup
    MEDDLY::unpacked_node::recycle(A);
    // Save in cache
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, value, result);
    return result;
 }

/////////////////////////////////////////////////////////////////////////////////////////

lv_selector_table::lv_selector_table(MEDDLY::expert_forest* forest)
: forest(forest)
{ }

lv_selector_op* lv_selector_table::get_op(size_t level) const
{
    if (table.empty())
        table.resize(forest->getNumVariables() + 1);
    if (table[level] == nullptr)
        table[level] = LV_SELECTOR_OPNAME->buildOperation(forest, forest, level);

    return table[level];
}

/////////////////////////////////////////////////////////////////////////////////////////



























/////////////////////////////////////////////////////////////////////////////////////////
// Enumerate all the distinct values that appear at a specified MDD level
/////////////////////////////////////////////////////////////////////////////////////////

void domain_values_enumerator::visit(const MEDDLY::expert_forest *argF,
                                     const MEDDLY::node_handle a) 
{
    if (a <= 0)
        return; // terminal

    if (a >= visited.size()) { // resize to adapt
        size_t nsz = visited.size();
        while (a >= nsz) nsz *= 4;
        visited.resize(nsz);
    }
    if (visited[a])
        return; // already visited

    assert(a < visited.size());
    visited[a] = true;

    MEDDLY::unpacked_node* A = argF->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node* A = MEDDLY::unpacked_node::newFromNode(argF, a, false);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // argF->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);

    const size_t a_size = get_node_size(A);
    const bool a_full = A->isFull();
    const int a_level = argF->getNodeLevel(a);

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        if (a_level == level) { // store the found values
            if (a_val < 0) {
                if (-a_val >= neg.size()) { // resize to adapt
                    size_t nsz = neg.size();
                    while (-a_val >= nsz) nsz *= 2;
                    neg.resize(nsz);
                }
                assert(-a_val < neg.size());
                neg[-a_val] = true;
            }
            else {
                if (a_val >= pos.size()) { // resize to adapt
                    size_t nsz = pos.size();
                    while (a_val >= nsz) nsz *= 2;
                    pos.resize(nsz);
                }
                assert(a_val < pos.size());
                pos[a_val] = true;
            }
        }
        else if (a_level > level) { // recursive visit
            visit(argF, A->d(i));
        }
    }

    // Cleanup
    MEDDLY::unpacked_node::recycle(A);
}

/////////////////////////////////////////////////////////////////////////////////////////

void domain_values_enumerator::get_values(std::vector<int>& out) const
{
    out.clear();
    for (ssize_t i=neg.size()-1; i>=0; i--) {
        if (neg[i])
            out.push_back(-i);
    }
    for (size_t i=0; i<pos.size(); i++) {
        if (pos[i])
            out.push_back(i);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////






















/////////////////////////////////////////////////////////////////////////////////////////
// Determine the smallest degree in a set, by level
/////////////////////////////////////////////////////////////////////////////////////////

smallest_degree_opname::smallest_degree_opname() 
/**/ : unary_opname("SmallestDegreeFinder") { }

smallest_degree_op* 
smallest_degree_opname::buildOperation(MEDDLY::forest* arg, 
                                       const smallest_degree_table* p_table, 
                                       const degree_type degtype, size_t lambda)
{
    if (0==arg) return 0;
    
    return new smallest_degree_op(this, (MEDDLY::expert_forest*)arg, p_table, degtype, lambda);
}

/////////////////////////////////////////////////////////////////////////////////////////

smallest_degree_op::smallest_degree_op(MEDDLY::unary_opname* oc, 
                                       MEDDLY::expert_forest* arg,
                                       const smallest_degree_table* p_table,
                                       const degree_type degtype,
                                       const size_t lambda) 
/**/ : base_NtoI(oc, arg), p_table(p_table), degtype(degtype), lambda(lambda)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

int smallest_degree_op::compute(MEDDLY::node_handle a)
{
    if (0 == a) 
        return -1;
    if (-1 == a)
        return 0;
    
    int res;
    MEDDLY::ct_entry_key* CTsrch;
    MEDDLY::ct_entry_key* key = findResult(a, res);
    if (nullptr==key)
        return res;

    // Read the node and compute the degrees
    int smallest_degree = -1;
    MEDDLY::unpacked_node* A = argF->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node* A = MEDDLY::unpacked_node::newFromNode(argF, a, false);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // argF->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);

    const size_t a_size = get_node_size(A);
    const size_t a_level = A->getLevel();
    const bool a_full = A->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        int down_degree = compute(A->d(i));

        int edge_degree = down_degree;
        if (p_table->pivot_order->is_below_lambda(lambda, a_level)) { 
            // level contributes to the degree
            edge_degree += get_degree_of(a_val, degtype);
        }

        if (smallest_degree < 0)
            smallest_degree = edge_degree;
        else
            smallest_degree = min(smallest_degree, edge_degree);
    }

    // Cleanup
    MEDDLY::unpacked_node::recycle(A);
    // Save in cache
    res = smallest_degree;
    saveResult(key, a, res);
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////

smallest_degree_table::smallest_degree_table(MEDDLY::expert_forest* forest,
                                             const variable_order *pivot_order) 
/**/ : forest(forest), pivot_order(pivot_order)
{ }

smallest_degree_op* smallest_degree_table::get_op(size_t level, const degree_type degtype) const 
{
    if (table.empty())
        table.resize(forest->getNumVariables() + 1);
    if (table[level].empty())
        table[level].resize(2);
    size_t dt = degtype==degree_type::BY_VALUE ? 1 : 0;
    if (nullptr == table[level][dt]) {
        table[level][dt] = SMALLEST_DEGREE_OPNAME->buildOperation
                            (forest, this, degtype, level);
    }
    return table[level][dt];
}

////////////////////////////////////////////////////////////////////////////////////////


















/////////////////////////////////////////////////////////////////////////////////////////
// Enumerate all degrees in a set, by level
/////////////////////////////////////////////////////////////////////////////////////////

degree_finder_opname::degree_finder_opname() 
/**/ : unary_opname("DegreeFinder") { }

degree_finder_op* 
degree_finder_opname::buildOperation(MEDDLY::forest* arg, 
                                     const degree_finder_table* p_table, 
                                     const degree_type degtype, size_t lambda)
{
    if (0==arg) return 0;
    
    return new degree_finder_op(this, (MEDDLY::expert_forest*)arg, p_table, degtype, lambda);
}

/////////////////////////////////////////////////////////////////////////////////////////

degree_finder_op::degree_finder_op(MEDDLY::unary_opname* oc, 
                                       MEDDLY::expert_forest* arg,
                                       const degree_finder_table* p_table,
                                       const degree_type degtype,
                                       const size_t lambda) 
/**/ : base_NtoI(oc, arg), p_table(p_table), degtype(degtype), lambda(lambda)
{ 
}

/////////////////////////////////////////////////////////////////////////////////////////

int degree_finder_op::compute(MEDDLY::node_handle a)
{
    if (0 == a) 
        return 0; // {}
    if (-1 == a)
        return 1; // {0}
    
    lut_key res;
    MEDDLY::ct_entry_key* CTsrch;
    MEDDLY::ct_entry_key* key = findResult(a, res);
    if (nullptr==key)
        return res;

    // Read the node and accumulate the set of all degrees
    std::vector<int> degrees;
    MEDDLY::unpacked_node* A = argF->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node* A = MEDDLY::unpacked_node::newFromNode(argF, a, false);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // argF->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);

    const size_t a_size = get_node_size(A);
    const size_t a_level = A->getLevel();
    const bool a_full = A->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));

        lut_key down_lk = compute(A->d(i));
        const std::vector<int>& down_degrees = p_table->look_up(down_lk);

        int edge_incr = 0;
        if (p_table->pivot_order->is_below_lambda(lambda, a_level)) { 
            // level contributes to the degree
            edge_incr = get_degree_of(a_val, degtype);
        }

        for (int d : down_degrees)
            degrees.push_back(add_exact(d, edge_incr));
    }

    // Sort and make the list of degrees unique
    std::sort(degrees.begin(), degrees.end(), [](int i, int j){ return i>j;});
    degrees.erase(std::unique(degrees.begin(), degrees.end()), degrees.end());

    // Cleanup
    MEDDLY::unpacked_node::recycle(A);
    // Save in cache
    res = p_table->lut.insert(std::move(degrees));
    saveResult(key, a, res);
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////

degree_finder_table::degree_finder_table(MEDDLY::expert_forest* forest,
                                         const variable_order *pivot_order) 
/**/ : forest(forest), pivot_order(pivot_order)
{
    lut.insert(std::vector<int>{});  // 0 terminal
    lut.insert(std::vector<int>{0});  // 1 terminal

    const size_t num_vars = forest->getNumVariables();
    table.resize(num_vars + 1);
    for (size_t lvl=0; lvl<=num_vars; ++lvl) {
        table[lvl].resize(2);
    }
}

degree_finder_op* degree_finder_table::get_op(size_t level, const degree_type degtype) const {
    size_t dt = degtype==degree_type::BY_VALUE ? 1 : 0;
    if (nullptr == table[level][dt]) {
        table[level][dt] = DEGREE_FINDER_OPNAME->buildOperation
                            (forest, this, degtype, level);
    }
    return table[level][dt];
}

////////////////////////////////////////////////////////////////////////////////////////



















/////////////////////////////////////////////////////////////////////////////////////////
// Select all elements that have a specific degree
/////////////////////////////////////////////////////////////////////////////////////////

degree_selector_opname::degree_selector_opname() 
/**/ : unary_opname("DegreeSelector")
{ }

degree_selector_op* 
degree_selector_opname::buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF,
                                       const degree_selector_table* p_table, 
                                       const degree_type degtype,
                                       size_t lambda)
{
    if (0==arF || 0==resF) return 0;
    
    return new degree_selector_op(this, (MEDDLY::expert_forest*)arF, 
                                  (MEDDLY::expert_forest*)resF,
                                  p_table, degtype, lambda);
}

degree_selector_op::~degree_selector_op() { }

/////////////////////////////////////////////////////////////////////////////////////////

degree_selector_op::degree_selector_op(MEDDLY::opname* opcode, MEDDLY::expert_forest* _argF,
                                       MEDDLY::expert_forest* _resF,
                                       const degree_selector_table* p_table,
                                       const degree_type degtype, const size_t lambda)
/**/ : base_NItoN(opcode, _argF, _resF), p_table(p_table), degtype(degtype), lambda(lambda)
{ }


/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::node_handle degree_selector_op::compute(MEDDLY::node_handle a, const int degree) 
{
    if (a==0 || a==-1) return (degree==0 ? a : 0);
    if (degree < 0) return 0;

    MEDDLY::node_handle result;
    MEDDLY::ct_entry_key* CTsrch;
    MEDDLY::ct_entry_key* key = findResult(a, degree, result);
    if (nullptr==key)
        return result;

    // // Read the node and accumulate the LCM of the GCDs
    MEDDLY::unpacked_node* A = argF->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    // MEDDLY::unpacked_node* A = MEDDLY::unpacked_node::newFromNode(argF, a, false);
    // MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New();
    // argF->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);
    const int a_level = argF->getNodeLevel(a);
    const size_t a_size = get_node_size(A);

    const int res_size = a_size;
    check_level_bound(resF, a_level, res_size);
    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, a_level, res_size);

    const bool a_full = A->isFull();

    for (size_t i = 0; i < (a_full ? a_size : A->getNNZs()); i++) { // for each a
        if (a_full && 0==A->d(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->i(i));
        
        int r_val = a_val;
        int down_degree = degree;
        if (p_table->pivot_order->is_below_lambda(lambda, a_level)) { 
            // level contributes to the degree
            down_degree -= get_degree_of(a_val, degtype);
        }

        unionNodes(C, compute(A->d(i), down_degree), 
                   ZtoNode(r_val), resF, mddUnion);
    }

    // // Cleanup
    MEDDLY::unpacked_node::recycle(A);
    // Save in cache
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, degree, result);
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

degree_selector_table::degree_selector_table(MEDDLY::expert_forest* forest,
                                             const variable_order *pivot_order) 
/**/ : forest(forest), pivot_order(pivot_order)
{
    const size_t num_vars = forest->getNumVariables();
    table.resize(num_vars + 1);
    for (size_t lvl=0; lvl<=num_vars; ++lvl) {
        table[lvl].resize(2);
    }
}

degree_selector_op* degree_selector_table::get_op(size_t level, const degree_type degtype) const
{
    size_t dt = degtype==degree_type::BY_VALUE ? 1 : 0;
    if (nullptr == table[level][dt]) {
        table[level][dt] = DEGREE_SELECTOR_OPNAME->buildOperation
                            (forest, forest, this, degtype, level);
    }
    return table[level][dt];
}

/////////////////////////////////////////////////////////////////////////////////////////





















/////////////////////////////////////////////////////////////////////////////////////////

void init_custom_meddly_operators(MEDDLY::forest* forestMDD, const variable_order *pivot_order) 
{
    MEDDLY::expert_forest *forestMDDexp = (MEDDLY::expert_forest*)forestMDD;

    // Declare out custom operators
    S_VECTORS_OPNAME = new s_vectors_opname();
    S_VECTORS_OPS = new s_vectors_table(forestMDDexp, pivot_order);

    LESSEQ_SQ_OPNAME = new lesseq_sq_opname();
    LESSEQ_SQ_OPS = new lesseq_sq_table(forestMDDexp, pivot_order);

    // COMPL_PROC_OPNAME = new compl_proc_opname();
    // COMPL_PROC_OPS = new compl_proc_table(forestMDDexp);

    SIGN_CANON_OPNAME = new sign_canon_mdd_opname();
    SIGN_CANON_OPS = new sign_canon_table(forestMDDexp);

    VMULT_OPNAME = new vmult_opname();
    VMULT = (vmult_op*)VMULT_OPNAME->buildOperation(forestMDDexp, forestMDD);

    VCANON_DIVISORS_MDD_OPNAME = new vcanon_mdd_opname(false);
    VCANON_DIVISORS_MDD = (vcanon_mdd_op*)VCANON_DIVISORS_MDD_OPNAME->buildOperation(forestMDD, forestMDD);

    VCANON_DIVIDE_MDD_OPNAME = new vcanon_mdd_opname(true);
    VCANON_DIVIDE_MDD = (vcanon_mdd_op*)VCANON_DIVIDE_MDD_OPNAME->buildOperation(forestMDD, forestMDD);

    DIV_FINDER_MDD_OPNAME = new divisors_finder_mdd_opname();
    DIV_FINDER_MDD = (divisors_finder_mdd_op*)DIV_FINDER_MDD_OPNAME->buildOperation(forestMDD);

    SUPPORT_INCL_OPNAME = new support_inclusion_opname();
    SUPPORT_INCL_TABLE = new support_inclusion_table(forestMDDexp, pivot_order);

    LV_SELECTOR_OPNAME = new lv_selector_opname();
    LV_SELECTOR_OPS = new lv_selector_table(forestMDDexp);

    SMALLEST_DEGREE_OPNAME = new smallest_degree_opname();
    SMALLEST_DEGREE_TABLE = new smallest_degree_table(forestMDDexp, pivot_order);

    DEGREE_FINDER_OPNAME = new degree_finder_opname();
    DEGREE_FINDER_TABLE = new degree_finder_table(forestMDDexp, pivot_order);

    DEGREE_SELECTOR_OPNAME = new degree_selector_opname();
    DEGREE_SELECTOR_TABLE = new degree_selector_table(forestMDDexp, pivot_order);
}

/////////////////////////////////////////////////////////////////////////////////////////

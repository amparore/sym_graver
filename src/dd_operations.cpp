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
            const char* var_name = mp.dd.getForest()->getDomain()->getVar(lvl + 1)->getName().c_str();
            col_spaces[var] = max(col_spaces[var], strlen(var_name));
            os << setw(col_spaces[var]) << var_name << " ";
        }
        os << endl;
    }
    // print paths/tuples/rows
    for (MEDDLY::enumerator path(mp.dd); path != 0; ++path) {
        const int *assignments = path.getAssignments() + 1;
        int degree = 0;
        for(size_t var=0; var < m; var++) {
            bool on = false;
            const size_t lvl = mp.vorder.var2lvl(var);
            int value = NodeToZ(assignments[lvl]);
            if (mp.lambda>0) {
                if (mp.pivot_order->is_above_lambda(mp.lambda, lvl+1)) { os << "\033[90m"; on = true; }
                else if (mp.lambda == lvl+1) { 
                    os << (value>0 ? "\033[92m" : (value==0 ? "\033[93m" : "\033[95m"));
                    on = true; 
                }
            }
            os << setw(col_spaces[var]) << right << value << " ";
            if (on) {
                os << "\033[39m";
            }
            if (mp.lambda>0 && mp.pivot_order->is_below_lambda(mp.lambda, lvl+1))
                degree += abs(value);
        }
        if (mp.lambda>0)
            os << " [" << degree << "]";
        os << endl;
    }
    return os;
}

/////////////////////////////////////////////////////////////////////////////////////////

// resize level to accomodate a given bound
void check_level_bound(MEDDLY::forest* forest, int level, int bound)
{
    // if (bound > 20)
    //     throw std::exception();

    if (forest->getDomain()->getVariableBound(level) < bound) {
        // cout << "resizing bound of level "<<level<<" from "<<dom_exp->getVariableBound(level)
        //      << " to "<<bound<<endl;
        forest->getDomain()->enlargeVariableBound(level, false, bound);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge 
mdd_from_vectors(const std::vector<std::vector<int>>& vecs, 
                 MEDDLY::forest* forestMDD, bool make_symmetric) 
{
    const size_t num_rows = make_symmetric ? 2*vecs.size() : vecs.size();
    MEDDLY::dd_edge set(forestMDD);
    if (num_rows) {
        std::vector<std::vector<int>> ins_data(num_rows);
        std::vector<const int*> ins_ptr;
        ins_ptr.resize(num_rows);

        for (size_t i=0; i<num_rows; i++) {
            size_t row = i % vecs.size();
            int sign = (i < vecs.size()) ? 1 : -1;
            assert(forestMDD->getDomain()->getNumVariables() == vecs.at(row).size());
            ins_data[i].resize(forestMDD->getDomain()->getNumVariables() + 1);

            // encode integers and verify bounds
            for (size_t j=0; j<vecs[row].size(); j++) {
                ins_data[i][j+1] = ZtoNode(sign * vecs[row][j]); // encode integers
                if (ins_data[i][j+1] >= forestMDD->getDomain()->getVariableBound(j+1)) {
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
    int bound = forestMDD->getDomain()->getVariableBound(level);
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
    int bound = forestMDD->getDomain()->getVariableBound(level);
    bool terms[bound+1];
    for (int i=0; i<bound+1; i++)
        terms[i] = (NodeToZ(i) == value);

    MEDDLY::dd_edge level_sel(forestMDD);
    forestMDD->createEdgeForVar(level, false, terms, level_sel);
    return level_sel;
}

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::dd_edge
selector_for_nonzeroes_at_level(MEDDLY::forest *forestMDD, size_t level) 
{
    int bound = forestMDD->getDomain()->getVariableBound(level);
    bool terms[bound+1];
    for (int i=0; i<bound+1; i++)
        terms[i] = (NodeToZ(i) != 0);

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
    MEDDLY::dd_edge zero(forestMDD);
    int num_levels = forestMDD->getDomain()->getNumVariables();
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
    assert(!n->isFull());
    if (n->isFull()) {
        // find exact sizes (exclude trailing zeroes)
        unsigned size = n->getSize();
        while (size > 0 && n->down(size - 1) == 0) 
            --size;
        return size;
    }
    else {
        // return the size containing the last nonzero
        return n->index( n->getSize()-1 ) + 1;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

void unionNodes(MEDDLY::unpacked_node* C, MEDDLY::node_handle node, 
                unsigned pos, MEDDLY::forest* resF,  
                MEDDLY::binary_operation* mddUnion) 
{
    if (node == 0)
        return; // nothing to add
    assert(pos < C->getSize());
    assert(C->isFull());
    if (C->down(pos) == 0)
        C->setFull(pos, node); // resF->linkNode(node); // FIXME: check
    else { // perform union
        MEDDLY::dd_edge dd1(resF), dd2(resF), union_dd(resF);
        dd1.set(node);
        dd2.set(C->down(pos));
        mddUnion->computeTemp(dd1, dd2, union_dd);
        assert(pos < resF->getDomain()->getVariableBound(C->getLevel()));
        C->setFull(pos, union_dd);
        assert(resF->getNodeLevel(C->down(pos)) <= C->getLevel());
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

void differenceNodes(MEDDLY::unpacked_node* C, MEDDLY::node_handle node, 
                     unsigned pos, MEDDLY::forest* resF,  
                     MEDDLY::binary_operation* mddDifference) 
{
    if (node == 0)
        return; // nothing to subtract
    assert(pos < C->getSize());
    assert(C->isFull());
    
    if (C->down(pos) != 0) { // perform difference
        MEDDLY::dd_edge dd1(resF), dd2(resF), diff_dd(resF);
        dd1.set(C->down(pos));
        dd2.set(node);
        mddDifference->computeTemp(dd1, dd2, diff_dd);
        assert(pos < resF->getDomain()->getVariableBound(C->getLevel()));
        C->setFull(pos, diff_dd);
        assert(resF->getNodeLevel(C->down(pos)) <= C->getLevel());
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

double dd_cardinality(const MEDDLY::dd_edge& dd) {
    double card;
    MEDDLY::apply(MEDDLY::CARDINALITY, dd, card);
    return card;
}

/////////////////////////////////////////////////////////////////////////////////////////














/////////////////////////////////////////////////////////////////////////////////////////

void dd_stats::generate_stats(const MEDDLY::dd_edge& e) {
    forest = e.getForest();
    const int max_levels = abs(e.getLevel());
    domain_values_per_lvl.resize(max_levels + 1);
    num_nodes_per_lvl.resize(max_levels + 1);
    num_edges_per_lvl.resize(max_levels + 1);
    visit(e.getNode());
    visited.clear();
}

/////////////////////////////////////////////////////////////////////////////////////////

void dd_stats::visit(const MEDDLY::node_handle node) {
    if (visited.count(node) > 0)
        return;
    if (node == forest->handleForValue(false) || 
        node == forest->handleForValue(true))
        return; // terminal node
    visited.insert(node);
    const int node_level = forest->getNodeLevel(node);

    MEDDLY::unpacked_node *rnode = MEDDLY::unpacked_node::New(forest);
    forest->unpackNode(rnode, node, MEDDLY::FULL_ONLY);
    assert(rnode->isFull());

    // get stats for this node
    ++num_nodes_per_lvl[node_level];
    for (int i = 0; i < rnode->getSize(); i++) {
        if (rnode->down(i) != forest->handleForValue(false)) {
            ++num_edges_per_lvl[node_level];
            int a_val = NodeToZ(i);

            if (domain_values_per_lvl[node_level].count(a_val) == 0)
                domain_values_per_lvl[node_level][a_val] = 0;

            ++domain_values_per_lvl[node_level][a_val];

            visit(rnode->down(i)); // recurse
        }
    }

    MEDDLY::unpacked_node::Recycle(rnode);
}

/////////////////////////////////////////////////////////////////////////////////////////

void dd_stats::write(std::ostream& os) const {
    // for (int i=1; i<num_nodes_per_lvl.size(); i++) {
    //     os << "Lvl " << i << ": ";
    //     os << "nodes: " << num_nodes_per_lvl[i] << " ";
    //     os << "edges: " << num_edges_per_lvl[i] << " ";
    //     os << "ratio: " << (double(num_edges_per_lvl[i]) / num_nodes_per_lvl[i]);

    //     os << " dom: [";
    //     for (auto v : domain_values_per_lvl[i]) {
    //         os << " " << v.first;
    //     }
    //     os << " ]";

    //     os << std::endl;
    // }
    os << (num_nodes_per_lvl.size() - 1) << std::endl;
    for (int i=1; i<num_nodes_per_lvl.size(); i++) {
        os << num_nodes_per_lvl[i] << " ";
        os << num_edges_per_lvl[i] << " ";
        os << domain_values_per_lvl[i].size() << " ";
        for (auto v : domain_values_per_lvl[i])
            os << v.first << " " << v.second << " ";
        os << std::endl;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

















/////////////////////////////////////////////////////////////////////////////////////////
// Base class for binary operators performing: Node * Node -> Node 
/////////////////////////////////////////////////////////////////////////////////////////

base_NNtoN::base_NNtoN(MEDDLY::binary_opname* opcode, MEDDLY::forest* arg1, 
                       MEDDLY::forest* arg2, MEDDLY::forest* res)
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

base_NItoN::base_NItoN(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
                       MEDDLY::forest* _resF)
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

base_NtoN::base_NtoN(MEDDLY::unary_opname* opcode, MEDDLY::forest* arg1, 
                     MEDDLY::forest* res)
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

base_NtoI::base_NtoI(MEDDLY::unary_opname* opcode, MEDDLY::forest* arg1)
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
// Base class for binary operators performing: Node * Node * Integer -> Node 
/////////////////////////////////////////////////////////////////////////////////////////

base_NNItoN::base_NNItoN(MEDDLY::opname* opcode, MEDDLY::forest* arg1, 
                         MEDDLY::forest* arg2, MEDDLY::forest* res)
  : MEDDLY::operation(opcode, 1), arg1F(arg1), arg2F(arg2), resF(res)
{
    MEDDLY::ct_entry_type* et;
    et = new MEDDLY::ct_entry_type(opcode->getName(), "NNI:N");
    et->setForestForSlot(0, arg1F);
    et->setForestForSlot(1, arg2F);
    et->setForestForSlot(4, resF);
    registerEntryType(0, et);
    buildCTs();

    mddUnion = MEDDLY::getOperation(MEDDLY::UNION, res, res, res);
}

base_NNItoN::~base_NNItoN() {
    unregisterInForest(arg1F);
    unregisterInForest(arg2F);   
    unregisterInForest(resF); 
}

bool base_NNItoN::checkForestCompatibility() const {
    return arg1F==resF;
}

/////////////////////////////////////////////////////////////////////////////////////////
 
inline MEDDLY::ct_entry_key* 
base_NNItoN::findResult(MEDDLY::node_handle a, MEDDLY::node_handle b, 
                        int i, MEDDLY::node_handle &c) 
{
    MEDDLY::ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    assert(CTsrch);
    // if (can_commute && a > b) {
    //     CTsrch->writeN(b);
    //     CTsrch->writeN(a);
    // } else {
        CTsrch->writeN(a);
        CTsrch->writeN(b);
    // }
    CTsrch->writeI(i);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    c = resF->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
}

inline void 
base_NNItoN::saveResult(MEDDLY::ct_entry_key* key, 
                        // MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
                        MEDDLY::node_handle c) 
{
    CTresult[0].reset();
    CTresult[0].writeN(c);
    CT0->addEntry(key, CTresult[0]);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Base of unary operations with signature: node * integer -> node * node
/////////////////////////////////////////////////////////////////////////////////////////

base_NItoNN::base_NItoNN(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
                         MEDDLY::forest* _res1F, MEDDLY::forest* _res2F)
/**/ : MEDDLY::operation(opcode, 1), argF(_argF), res1F(_res1F), res2F(_res2F)
{
    registerInForest(argF);
    registerInForest(res1F);
    registerInForest(res2F);

    MEDDLY::ct_entry_type* et;
    et = new MEDDLY::ct_entry_type(opcode->getName(), "NI:NN");
    et->setForestForSlot(0, argF);
    et->setForestForSlot(3, res1F);
    et->setForestForSlot(4, res2F);
    registerEntryType(0, et);
    buildCTs();

    mddUnion = MEDDLY::getOperation(MEDDLY::UNION, res1F, res1F, res1F);  
}

base_NItoNN::~base_NItoNN() {
    unregisterInForest(argF);
    unregisterInForest(res1F);   
    unregisterInForest(res2F);   
}

////////////////////////////////////////////////////////////////////////////////////////

inline void 
base_NItoNN::saveResult(MEDDLY::ct_entry_key* key, 
                        //MEDDLY::node_handle a, const int b, 
                        std::pair<MEDDLY::node_handle, MEDDLY::node_handle> c) 
{
    CTresult[0].reset();
    CTresult[0].writeN(c.first);
    CTresult[0].writeN(c.second);
    CT0->addEntry(key, CTresult[0]);
}

/////////////////////////////////////////////////////////////////////////////////////////

inline MEDDLY::ct_entry_key* 
base_NItoNN::findResult(MEDDLY::node_handle a, const int b, 
                        std::pair<MEDDLY::node_handle, MEDDLY::node_handle> &c) 
{
    MEDDLY::ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    assert(CTsrch);
    CTsrch->writeN(a);
    CTsrch->writeI(b);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    c.first  = res1F->linkNode(CTresult[0].readN());
    c.second = res2F->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////

bool base_NItoNN::checkForestCompatibility() const {
    return argF==res1F;
}

/////////////////////////////////////////////////////////////////////////////////////////

void base_NItoNN::computeDDEdge(const MEDDLY::dd_edge &a, const int b, 
                                MEDDLY::dd_edge &res1, MEDDLY::dd_edge &res2)
{
    std::pair<MEDDLY::node_handle, MEDDLY::node_handle> cnodes;
    cnodes = compute(a.getNode(), b);
//   const int num_levels = resF->getDomain()->getNumVariables();
//   if ( userFlag && resF->isQuasiReduced() && cnode != resF->getTransparentNode()
//     && resF->getNodeLevel(cnode) < num_levels) {
//     node_handle temp = ((mt_forest*)resF)->makeNodeAtLevel(num_levels, cnode);
//     resF->unlinkNode(cnode);
//     cnode = temp;
//   }
    res1.set(cnodes.first);
    res2.set(cnodes.second);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Base class for binary operators performing: Node * Node * Integer -> Node * Node
/////////////////////////////////////////////////////////////////////////////////////////

base_NNItoNN::base_NNItoNN(MEDDLY::opname* opcode, 
                           MEDDLY::forest* arg1, MEDDLY::forest* arg2, 
                           MEDDLY::forest* res1, MEDDLY::forest* res2)
: MEDDLY::operation(opcode, 1), arg1F(arg1), arg2F(arg2), res1F(res1), res2F(res2)
{
    MEDDLY::ct_entry_type* et;
    et = new MEDDLY::ct_entry_type(opcode->getName(), "NNI:NN");
    et->setForestForSlot(0, arg1F);
    et->setForestForSlot(1, arg2F);
    et->setForestForSlot(4, res1F);
    et->setForestForSlot(5, res2F);
    registerEntryType(0, et);
    buildCTs();

    mddUnion = MEDDLY::getOperation(MEDDLY::UNION, arg1F, arg2F, res1F);
}

base_NNItoNN::~base_NNItoNN() {
    unregisterInForest(arg1F);
    unregisterInForest(arg2F);   
    unregisterInForest(res1F); 
    unregisterInForest(res2F); 
}

bool base_NNItoNN::checkForestCompatibility() const {
    return arg1F==res1F;
}

/////////////////////////////////////////////////////////////////////////////////////////
 
inline MEDDLY::ct_entry_key* 
base_NNItoNN::findResult(MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
                         std::pair<MEDDLY::node_handle, MEDDLY::node_handle> &c) 
{
    MEDDLY::ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    assert(CTsrch);
    CTsrch->writeN(a);
    CTsrch->writeN(b);
    CTsrch->writeI(i);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    c.first  = res1F->linkNode(CTresult[0].readN());
    c.second = res2F->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
}

inline void 
base_NNItoNN::saveResult(MEDDLY::ct_entry_key* key, 
                         //MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
                         std::pair<MEDDLY::node_handle, MEDDLY::node_handle> c)
{
    CTresult[0].reset();
    CTresult[0].writeN(c.first);
    CTresult[0].writeN(c.second);
    CT0->addEntry(key, CTresult[0]);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Base class for binary operators performing: Node * Node * Integer -> Node * Node * Node
/////////////////////////////////////////////////////////////////////////////////////////

base_NNItoNNN::base_NNItoNNN(MEDDLY::opname* opcode, 
                             MEDDLY::forest* arg1, MEDDLY::forest* arg2, 
                             MEDDLY::forest* res1, MEDDLY::forest* res2,
                             MEDDLY::forest* res3)
: MEDDLY::operation(opcode, 1), arg1F(arg1), arg2F(arg2), 
  res1F(res1), res2F(res2), res3F(res3)
{
    MEDDLY::ct_entry_type* et;
    et = new MEDDLY::ct_entry_type(opcode->getName(), "NNI:NNN");
    et->setForestForSlot(0, arg1F);
    et->setForestForSlot(1, arg2F);
    et->setForestForSlot(4, res1F);
    et->setForestForSlot(5, res2F);
    et->setForestForSlot(6, res3F);
    registerEntryType(0, et);
    buildCTs();
}

base_NNItoNNN::~base_NNItoNNN() {
    unregisterInForest(arg1F);
    unregisterInForest(arg2F);   
    unregisterInForest(res1F); 
    unregisterInForest(res2F); 
    unregisterInForest(res3F); 
}

bool base_NNItoNNN::checkForestCompatibility() const {
    return arg1F==res1F;
}

/////////////////////////////////////////////////////////////////////////////////////////
 
inline MEDDLY::ct_entry_key* 
base_NNItoNNN::findResult(MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
                          std::tuple<MEDDLY::node_handle, MEDDLY::node_handle, MEDDLY::node_handle> &c) 
{
    MEDDLY::ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
    assert(CTsrch);
    CTsrch->writeN(a);
    CTsrch->writeN(b);
    CTsrch->writeI(i);
    CT0->find(CTsrch, CTresult[0]);
    if (!CTresult[0]) return CTsrch;
    std::get<0>(c) = res1F->linkNode(CTresult[0].readN());
    std::get<1>(c) = res2F->linkNode(CTresult[0].readN());
    std::get<2>(c) = res3F->linkNode(CTresult[0].readN());
    CT0->recycle(CTsrch);
    return 0;
}

inline void 
base_NNItoNNN::saveResult(MEDDLY::ct_entry_key* key, 
                         //MEDDLY::node_handle a, MEDDLY::node_handle b, int i, 
                         std::tuple<MEDDLY::node_handle, MEDDLY::node_handle, MEDDLY::node_handle> c)
{
    CTresult[0].reset();
    CTresult[0].writeN(std::get<0>(c));
    CTresult[0].writeN(std::get<1>(c));
    CTresult[0].writeN(std::get<2>(c));
    CT0->addEntry(key, CTresult[0]);
}

/////////////////////////////////////////////////////////////////////////////////////////













/////////////////////////////////////////////////////////////////////////////////////////
// S-Vectors
/////////////////////////////////////////////////////////////////////////////////////////

typedef union s_vector_flags_t { // cache entry
    int value;
    struct {
        // the two summed vectors are potentially conformant, i.e. there are no
        // opposite signs in any entries of the vectors. A sum of conformant vectors is
        // surely recudible, and does not need to be generated
        bool is_potentially_conformant : 1; 
        // Compute (a-b) instead of (a+b)
        bool subtract_b : 1;
        // When computing the Graver half-basis, we need a canonical form for the vectors
        // We decide arbitrarily that the first entry from the DD root has to be positive.
        // If it is negative, the resulting vector (a+b) is inverted as (-(a+b))
        unsigned sign_of_sum : 2;
        // The level (for project-and-lift)
        unsigned lambda : (sizeof(int)*8 - 4);
    } bf;
} s_vector_flags_t;

static_assert(sizeof(s_vector_flags_t) == sizeof(int));

/////////////////////////////////////////////////////////////////////////////////////////

s_vectors::s_vectors(MEDDLY::opname* opcode, MEDDLY::forest* arg1, 
                     MEDDLY::forest* arg2, MEDDLY::forest* res,
                     const variable_order *pivot_order)
: base_NNItoN(opcode, arg1, arg2, res), pivot_order(pivot_order)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::node_handle s_vectors::compute(MEDDLY::node_handle a, MEDDLY::node_handle b, 
                                       const bool is_potentially_conformant, 
                                       const ab_sum_t sum_or_diff, 
                                       sv_sign sign_of_sum, const size_t lambda) 
{
    if (a==0 || b==0) return 0;
    if (a==-1 && b==-1) return is_potentially_conformant ? 0 : -1;

    // Commutativity (only for sum)
    if (sum_or_diff==ab_sum_t::A_PLUS_B && a < b)
        std::swap(a, b);

    if (sum_or_diff==ab_sum_t::A_MINUS_B && sign_of_sum!=sv_sign::SVS_UNDECIDED && a < b) {
        // cout << "*" <<endl;
        std::swap(a, b);
        sign_of_sum = (sign_of_sum==sv_sign::SVS_POS) ? sv_sign::SVS_NEG : sv_sign::SVS_POS;
    }

    s_vector_flags_t svf;
    svf.bf.is_potentially_conformant = is_potentially_conformant;
    svf.bf.subtract_b = (sum_or_diff==ab_sum_t::A_MINUS_B ? 1 : 0);
    svf.bf.sign_of_sum = sign_of_sum;
    svf.bf.lambda = lambda;

    MEDDLY::node_handle result;
    MEDDLY::ct_entry_key* key = findResult(a, b, svf.value, result);
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

    const size_t a_size = get_node_size(A);
    const size_t b_size = get_node_size(B);
    const int res_size = a_size + b_size + 1; // with +1 because we are encoding integers
    check_level_bound(resF, res_level, res_size);
    
    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, res_level, res_size);

    const bool a_full = A->isFull(), b_full = B->isFull();

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        for (unsigned j = 0; j < (b_full ? b_size : B->getSize()); j++) { // for each b
            if (b_full && 0==B->down(j))
                continue;
            int b_val = NodeToZ(b_full ? j : B->index(j));

            if (sum_or_diff==ab_sum_t::A_MINUS_B)
                b_val = -b_val;

            int ab_sum = add_exact(a_val, b_val); // a + b or a - b
            // int ab_sign_prod = sign3(a_val) * sign3(b_val); // a * b
            // canonicalize the sum (if not yet decided)
            sv_sign curr_sign_of_sum = sign_of_sum;
            if (sign_of_sum == SVS_UNDECIDED && 0!=ab_sum) { 
                curr_sign_of_sum = (ab_sum > 0) ? SVS_POS : SVS_NEG;
            }
            // compute -(a+b) or -(a-b)
            if (curr_sign_of_sum == SVS_NEG) {
                ab_sum = -ab_sum;
            }

            int ab_sum_idx = ZtoNode(ab_sum);
            if (ab_sum_idx >= res_size)
                cout << "ab_sum_idx:"<<ab_sum_idx<<" res_size:"<<res_size<<endl;
            assert(ab_sum_idx < res_size);

            bool ij_conf = comparable_signs(a_val, b_val);

            bool down_is_conf = is_potentially_conformant && ij_conf;

            bool do_sum = true;
            if (lambda==0) {
            }
            else if (pivot_order->is_same_as_lambda(lambda, res_level)) {
                do_sum = !ij_conf;
            }
            else if (pivot_order->is_below_lambda(lambda, res_level)) {
                do_sum = ij_conf;
            }
            else if (pivot_order->is_above_lambda(lambda, res_level)) {
            }
            // cout << "a_level:"<<a_level<<" a_val:"<<a_val<<" b_val:"<<b_val<<" ij_conf:"<<ij_conf<<" do_sum:"<<do_sum<<endl;

            if (do_sum) {
                MEDDLY::node_handle sum;
                sum = compute(A->down(i), B->down(j), down_is_conf, sum_or_diff, curr_sign_of_sum, lambda);
                unionNodes(C, sum, ab_sum_idx, resF, mddUnion);
            }
        }
    }

    // cleanup
    MEDDLY::unpacked_node::Recycle(B);
    MEDDLY::unpacked_node::Recycle(A);

    // reduce and return result
    result = resF->createReducedNode(-1, C);
    saveResult(key, /*a, b, svf.value,*/ result);
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

void 
s_vectors::computeDDEdge(const MEDDLY::dd_edge &ar1, const MEDDLY::dd_edge &ar2, 
                         const bool is_potentially_conformant, const ab_sum_t sum_or_diff, 
                         const sv_sign sign_of_sum, const size_t lambda,
                         MEDDLY::dd_edge &res)
{
    MEDDLY::node_handle cnode = compute(ar1.getNode(), ar2.getNode(),
                            is_potentially_conformant, sum_or_diff, sign_of_sum, lambda);
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

s_vectors* 
s_vectors_opname::buildOperation(MEDDLY::forest* a1, 
                                 MEDDLY::forest* a2, MEDDLY::forest* r,
                                 const variable_order *pivot_order)
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
        return new s_vectors(this, a1, a2, r, pivot_order);
    }
    throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

/////////////////////////////////////////////////////////////////////////////////////////












// /////////////////////////////////////////////////////////////////////////////////////////
// // Less Equal Squared Operator
// // Get the elements using the less-equal-but-not-equal-squared operator
// /////////////////////////////////////////////////////////////////////////////////////////

// lesseq_sq::lesseq_sq(MEDDLY::binary_opname* opcode, 
// /**/        MEDDLY::forest* arg1, MEDDLY::forest* arg2, MEDDLY::forest* res,
// /**/        const lesseq_sq_table* tab,
// /**/        bool isEq, bool isB0,
// /**/        bool _subtract, size_t lambda)
//   : base_NNtoN(opcode, arg1, arg2, res), p_table(tab),
//     isPotentiallyEqual(isEq), isBZero(isB0),
//     subtract(_subtract), lambda(lambda)
// { }

// /////////////////////////////////////////////////////////////////////////////////////////

// MEDDLY::node_handle lesseq_sq::compute(MEDDLY::node_handle a, MEDDLY::node_handle b) {
//     if (a==0) return 0;
//     if (b==0) return 0;
//     if (a==-1) return isPotentiallyEqual || isBZero ? 0 : -1;
//     assert(b != -1);

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
//         : arg1F->newUnpacked(a, MEDDLY::SPARSE_ONLY);
//     MEDDLY::unpacked_node *B = (b_level < res_level)
//         ? MEDDLY::unpacked_node::newRedundant(arg2F, res_level, b, false)
//         : arg2F->newUnpacked(b, MEDDLY::SPARSE_ONLY);

//     const size_t a_size = get_node_size(A);
//     const size_t b_size = get_node_size(B);
//     // const size_t res_size = a_size;
//     size_t res_size;
//     if (lambda != 0 && 
//         p_table->pivot_order->is_above_lambda(lambda, res_level) && subtract) 
//     {
//         res_size = a_size + b_size;
//     }
//     else res_size = a_size;
//     check_level_bound(resF, res_level, res_size);

//     MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, res_level, res_size);

//     const bool a_full = A->isFull(), b_full = B->isFull();

//     for (size_t i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
//         if (a_full && 0==A->down(i))
//             continue;
//         int a_val = NodeToZ(a_full ? i : A->index(i));

//         for (size_t j = 0; j < (b_full ? b_size : B->getSize()); j++) { // for each b
//             if (b_full && 0==B->down(j))
//                 continue;
//             int b_val = NodeToZ(b_full ? j : B->index(j));

//             bool is_lesseq_sq;
//             bool is_pot_eq = isPotentiallyEqual;
//             bool is_b_pot_zero = isBZero;
//             if (lambda != 0 && p_table->pivot_order->is_above_lambda(lambda, res_level)) {
//                 is_lesseq_sq = true;
//             }
//             else {
//             //     // check that i <= j and both are conformal
//                 int ab_sign_prod = sign3(a_val) * sign3(b_val);
//                 is_lesseq_sq = (abs(b_val) <= abs(a_val) && ab_sign_prod >= 0);
//                 is_pot_eq &= (a_val == b_val);
//                 is_b_pot_zero &= (0 == b_val);
//             }

//             if (is_lesseq_sq) 
//             { 
//                 int idx = ZtoNode(subtract ? subtract_exact(a_val, b_val) : a_val);

//                 // Determine next operator in chain
//                 bool down_is_eq = isPotentiallyEqual && is_pot_eq;
//                 bool down_is_b0 = isBZero && is_b_pot_zero;
//                 lesseq_sq* next_op = p_table->get_op(lambda, down_is_eq, 
//                                                      down_is_b0, subtract);

//                 MEDDLY::node_handle leq_ij = next_op->compute(A->down(i), B->down(j));
//                 unionNodes(C, leq_ij, idx, resF, mddUnion);
//             }
//         }
//     }

//     // cleanup
//     MEDDLY::unpacked_node::Recycle(B);
//     MEDDLY::unpacked_node::Recycle(A);
//     // reduce and return result
//     result = resF->createReducedNode(-1, C);
//     saveResult(key, a, b, result);
//     return result;
// }

// /////////////////////////////////////////////////////////////////////////////////////////

// lesseq_sq* 
// lesseq_sq_opname::buildOperation(MEDDLY::forest* a1, 
//                                           MEDDLY::forest* a2, 
//                                           MEDDLY::forest* r,
//                                           const lesseq_sq_table* tab,
//                                           bool isEq, bool isB0, bool subtract,
//                                           size_t lambda)
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
//         return new lesseq_sq(this, a1, a2, r, tab, isEq, isB0,
//                              subtract, lambda);
//     }
//     throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED, __FILE__, __LINE__);
// }

// MEDDLY::binary_operation* 
// lesseq_sq_opname::buildOperation(MEDDLY::forest* a1, 
//                                           MEDDLY::forest* a2, 
//                                           MEDDLY::forest* r)
// {
//     throw std::exception(); // Unimplemented
// }

// /////////////////////////////////////////////////////////////////////////////////////////

// lesseq_sq_table::lesseq_sq_table(MEDDLY::forest* forest,
//                                  const variable_order *pivot_order) 
// /**/ : forest(forest), pivot_order(pivot_order)
// { }

// /////////////////////////////////////////////////////////////////////////////////////////

// lesseq_sq* 
// lesseq_sq_table:: get_op(size_t level, bool isPotentiallyEqual, 
//                          bool isBZero, bool subtract) const 
// {
//     if (table.empty())
//         table.resize(forest->getNumVariables() + 1);
//     if (table[level].empty())
//         table[level].resize(2);
//     if (table[level][isPotentiallyEqual].empty())
//         table[level][isPotentiallyEqual].resize(2);
//     if (table[level][isPotentiallyEqual][isBZero].empty())
//         table[level][isPotentiallyEqual][isBZero].resize(2);
//     if (table[level][isPotentiallyEqual][isBZero][subtract] == nullptr)
//         table[level][isPotentiallyEqual][isBZero][subtract] = 
//             LESSEQ_SQ_OPNAME->buildOperation(forest, forest, forest, this, 
//                                              isPotentiallyEqual, isBZero, 
//                                              subtract, level);

//     return table[level][isPotentiallyEqual][isBZero][subtract];
// }

// /////////////////////////////////////////////////////////////////////////////////////////
















/////////////////////////////////////////////////////////////////////////////////////////
// Less-Equal-but-Not-Equal-Squared Operator
// Get the elements using the less-equal-but-not-equal-squared operator
/////////////////////////////////////////////////////////////////////////////////////////

leq_neq_sq::leq_neq_sq(MEDDLY::opname* opcode, MEDDLY::forest* forestMDD,
                       const variable_order *pivot_order, const bool subtract)
: base_NNItoN(opcode, forestMDD, forestMDD, forestMDD), 
/**/ pivot_order(pivot_order), subtract(subtract)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

typedef union leq_neq_sq_flags_t { // cache entry
    int value;
    struct {
        // is a and b potentially the same vectors?
        bool is_potentially_equal : 1; 
        // is b potentially the zero vector?
        bool is_b_potentially_zero : 1;
        // The level (for project-and-lift)
        unsigned lambda : (sizeof(int)*8 - 2);
    } bf;
} leq_neq_sq_flags_t;

static_assert(sizeof(leq_neq_sq_flags_t) == sizeof(int));

/////////////////////////////////////////////////////////////////////////////////////////

MEDDLY::node_handle 
leq_neq_sq::compute(MEDDLY::node_handle a, MEDDLY::node_handle b, int flags) 
{
    leq_neq_sq_flags_t lf { .value=flags };

    if (a==0) return 0;
    if (b==0) return 0;
    if (a==-1) return lf.bf.is_potentially_equal || lf.bf.is_b_potentially_zero ? 0 : -1;
    assert(b != -1);

    MEDDLY::node_handle result;
    MEDDLY::ct_entry_key* key = findResult(a, b, flags, result);
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

    const size_t a_size = get_node_size(A);
    const size_t b_size = get_node_size(B);
    // const size_t res_size = a_size;
    size_t res_size;
    if (subtract && lf.bf.lambda != 0 && 
        pivot_order->is_above_lambda(lf.bf.lambda, res_level)) 
    {
        res_size = a_size + b_size;
    }
    else res_size = a_size;
    check_level_bound(resF, res_level, res_size);

    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, res_level, res_size);

    const bool a_full = A->isFull(), b_full = B->isFull();

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        for (unsigned j = 0; j < (b_full ? b_size : B->getSize()); j++) { // for each b
            if (b_full && 0==B->down(j))
                continue;
            int b_val = NodeToZ(b_full ? j : B->index(j));

            bool ij_is_leq_sq;
            bool ij_is_pot_eq = lf.bf.is_potentially_equal;
            bool ij_is_b_pot_zero = lf.bf.is_b_potentially_zero;

            if (lf.bf.lambda != 0 && pivot_order->is_above_lambda(lf.bf.lambda, res_level)) {
                ij_is_leq_sq = true;
            }
            else {
            //     // check that i <= j and have comparable signs
                // int ab_sign_prod = sign3(a_val) * sign3(b_val);
                // ij_is_leq_sq = (abs(b_val) <= abs(a_val) && ab_sign_prod >= 0);
                ij_is_leq_sq = less_equal_squared(b_val, a_val);
                ij_is_pot_eq = ij_is_pot_eq && (a_val == b_val);
                ij_is_b_pot_zero = ij_is_b_pot_zero && (0 == b_val);
            }

            if (ij_is_leq_sq) 
            { 
                int idx = ZtoNode(subtract ? subtract_exact(a_val, b_val) : a_val);

                leq_neq_sq_flags_t down_lf;
                down_lf.bf.is_potentially_equal = ij_is_pot_eq;
                down_lf.bf.is_b_potentially_zero = ij_is_b_pot_zero;
                down_lf.bf.lambda = lf.bf.lambda;

                MEDDLY::node_handle leq_ij = compute(A->down(i), B->down(j), down_lf.value);
                unionNodes(C, leq_ij, idx, resF, mddUnion);
            }
        }
    }

    // cleanup
    MEDDLY::unpacked_node::Recycle(B);
    MEDDLY::unpacked_node::Recycle(A);
    // reduce and return result
    result = resF->createReducedNode(-1, C);
    saveResult(key, result);
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

void 
leq_neq_sq:: computeDDEdge(const MEDDLY::dd_edge &a, const MEDDLY::dd_edge &b, 
                           const bool is_potentially_equal, 
                           const bool is_b_potentially_zero,
                           const size_t lambda,
                           MEDDLY::dd_edge &res)
{
    leq_neq_sq_flags_t lf;
    lf.bf.is_potentially_equal = is_potentially_equal;
    lf.bf.is_b_potentially_zero = is_b_potentially_zero;
    lf.bf.lambda = lambda;

    MEDDLY::node_handle cnode = compute(a.getNode(), b.getNode(), lf.value);

    res.set(cnode);
}

/////////////////////////////////////////////////////////////////////////////////////////

leq_neq_sq* 
leq_neq_sq_opname::buildOperation(MEDDLY::forest* forestMDD,
                                  const variable_order *pivot_order,
                                  const bool subtract)
{
    if (0==forestMDD) return 0;

    if (forestMDD->isForRelations())
        throw MEDDLY::error(MEDDLY::error::TYPE_MISMATCH, __FILE__, __LINE__);

    if (forestMDD->getEdgeLabeling() == MEDDLY::edge_labeling::MULTI_TERMINAL) {
        if (forestMDD->isForRelations())
            throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED);
        return new leq_neq_sq(this, forestMDD, pivot_order, subtract);
    }
    throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

/////////////////////////////////////////////////////////////////////////////////////////









/////////////////////////////////////////////////////////////////////////////////////////
// Reduction of elements that are less-equal-squared-but-not-equal up to lambda
/////////////////////////////////////////////////////////////////////////////////////////

typedef union reduce_flags_t { // cache entry
    int value;
    struct {
        // is a and b potentially the same vectors?
        bool is_potentially_equal : 1; 
        // is b potentially the zero vector?
        bool is_b_potentially_zero : 1;
        // When computing the Graver half-basis, we need a canonical form for the redisual 
        // vectors. We decide arbitrarily that the first entry from the DD root has to be 
        // positive. If it is negative, the resulting vector (a-b) is inverted as (-(a-b))
        unsigned sign_of_sum : 2;
        // When comparing the Graver half-basis, we check for both b < a and -b < a,
        // and decide the comparison sign at the first occasion possible
        unsigned sign_of_comparison : 2;
        // The level (for project-and-lift)
        unsigned lambda : (sizeof(int)*8 - 6);
    } bf;
} reduce_flags_t;

static_assert(sizeof(reduce_flags_t) == sizeof(int));

/////////////////////////////////////////////////////////////////////////////////////////

reduce::reduce(MEDDLY::opname* opcode, MEDDLY::forest* forestMDD,
               const variable_order *pivot_order, const bool compute_differences)
: base_NNItoNN(opcode, forestMDD, forestMDD, forestMDD, forestMDD), 
  pivot_order(pivot_order), compute_differences(compute_differences)
{ 
    mddUnion = MEDDLY::getOperation(MEDDLY::UNION, forestMDD, forestMDD, forestMDD);
}

/////////////////////////////////////////////////////////////////////////////////////////

void 
reduce::computeDDEdge(const MEDDLY::dd_edge &a, const MEDDLY::dd_edge &b, 
                      const bool is_potentially_equal, 
                      const bool is_b_potentially_zero,
                      const sv_sign sign_of_sum, 
                      const cmp_sign sign_of_comparison,
                      const size_t lambda,
                      MEDDLY::dd_edge &irreducibles, 
                      MEDDLY::dd_edge &reduced) 
{
    reduce_flags_t rf;
    rf.bf.is_potentially_equal = is_potentially_equal;
    rf.bf.is_b_potentially_zero = is_b_potentially_zero;
    rf.bf.sign_of_sum = sign_of_sum;
    rf.bf.sign_of_comparison = sign_of_comparison;
    rf.bf.lambda = lambda;

    auto cnodes = compute(a.getNode(), b.getNode(), rf.value);

    irreducibles.set(cnodes.first);
    reduced.set(cnodes.second);
}

/////////////////////////////////////////////////////////////////////////////////////////

// first=irreducibles
// second=reduced
std::pair<MEDDLY::node_handle, MEDDLY::node_handle>
reduce::compute(MEDDLY::node_handle a, MEDDLY::node_handle b, const int flags) 
{
    reduce_flags_t rf;
    rf.value = flags;
    // counter_steps++;

    if (a==0) return std::make_pair(0, 0);
    if (b==0) return std::make_pair(res1F->linkNode(a), 0);
    if (a==-1) {
        if (rf.bf.is_potentially_equal || rf.bf.is_b_potentially_zero)
            return std::make_pair(-1, 0);
        else
            return std::make_pair(0, -1);
    }
    assert(b != -1);

    std::pair<MEDDLY::node_handle, MEDDLY::node_handle> result;
    MEDDLY::ct_entry_key* key = findResult(a, b, flags, result);
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

    const size_t a_size = get_node_size(A);
    const size_t b_size = get_node_size(B);
    const size_t resA_size = a_size;
    const size_t resAB_size = a_size + b_size + 1; // TODO: since b<=a, then a-b is <= a and reY_size should be a_size
    // if (rf.bf.lambda != 0 && 
    //     pivot_order->is_above_lambda(rf.bf.lambda, res_level)) 
    // {
    //     resAB_size = a_size + b_size;
    // }
    // else resAB_size = a_size;
    check_level_bound(res1F, res_level, resAB_size);
    MEDDLY::unpacked_node* C_irreducibles1;
    MEDDLY::unpacked_node* C_reduced2;

    C_irreducibles1 = MEDDLY::unpacked_node::newFull(res1F, res_level, resA_size);
    if (compute_differences)
        C_reduced2 = MEDDLY::unpacked_node::newFull(res2F, res_level, resAB_size);

    const bool a_full = A->isFull(), b_full = B->isFull();

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        MEDDLY::dd_edge irred_Ai(res1F);
        irred_Ai.set(res1F->linkNode(A->down(i)));

        for (unsigned j = 0; j < (b_full ? b_size : B->getSize()); j++) { // for each b
            if (b_full && 0==B->down(j))
                continue;
            int b_val = NodeToZ(b_full ? j : B->index(j));
            // if (abs(b_val) > abs(a_val))
            //     continue; // there are no more j that can reduce i.

            bool ij_reduce;
            bool ij_pot_eq = rf.bf.is_potentially_equal;
            bool ij_b_pot_zero = rf.bf.is_b_potentially_zero;
            cmp_sign ij_cmp_sign = (cmp_sign)rf.bf.sign_of_comparison;
            if (rf.bf.lambda != 0 && pivot_order->is_above_lambda(rf.bf.lambda, res_level)) {
                ij_reduce = true;
            }
            else {
                if (abs(b_val) > abs(a_val))
                    continue; // there are no more j that can reduce i.

                if (ij_cmp_sign == CMP_UNDECIDED) {
                    if      (b_val * a_val > 0)   ij_cmp_sign = CMP_POS;
                    else if (b_val * a_val < 0)   ij_cmp_sign = CMP_NEG;
                }
                int b_cmp_sign = (ij_cmp_sign == CMP_NEG ? -1 : 1);
                ij_reduce = less_equal_squared(b_val * b_cmp_sign, a_val);
                ij_pot_eq = ij_pot_eq && (a_val == b_val);
                ij_b_pot_zero = ij_b_pot_zero && (0 == b_val);
            }

            if (ij_reduce)
            {
                int b_cmp_sign = (ij_cmp_sign == CMP_NEG ? -1 : 1);
                // assert(b_val * a_val==0 || ij_cmp_sign!=CMP_UNDECIDED); // TODO: enable again
                int a_minus_b = subtract_exact(a_val, b_val * b_cmp_sign);
                // decide the sign of the result (if still undecided)
                sv_sign curr_sign_of_sum = (sv_sign)rf.bf.sign_of_sum;
                if (curr_sign_of_sum == SVS_UNDECIDED && 0!=a_minus_b) { 
                    curr_sign_of_sum = (a_minus_b > 0) ? SVS_POS : SVS_NEG;
                }
                // compute (a-b) or -(a-b)
                if (curr_sign_of_sum == SVS_NEG) {
                    a_minus_b = -a_minus_b;
                }

                reduce_flags_t down_rf;
                down_rf.bf.is_potentially_equal = ij_pot_eq;
                down_rf.bf.is_b_potentially_zero = ij_b_pot_zero;
                down_rf.bf.sign_of_sum = curr_sign_of_sum;
                down_rf.bf.sign_of_comparison = ij_cmp_sign;
                down_rf.bf.lambda = rf.bf.lambda;

                std::pair<MEDDLY::node_handle, MEDDLY::node_handle> down;
                down = compute(irred_Ai.getNode(), B->down(j), down_rf.value);

                // cout << "lvl(a)="<<a_level<<" lvl(b)="<<b_level<<" a="<<a_val<<" b="<<b_val
                //      << " is_eq="<<ij_pot_eq<<" b_zero="<<ij_b_pot_zero
                //      << " sign_of_sum="<<down_rf.bf.sign_of_sum
                //      << " sign_of_comparison="<<down_rf.bf.sign_of_comparison<<endl;

                irred_Ai.set(down.first);
                if (compute_differences)
                    unionNodes(C_reduced2, down.second, ZtoNode(a_minus_b), res2F, mddUnion);

                // THIS ONLY MAKES SENSE WHEN ALL PIVOT LEVELS ARE TO THE BOTTOM OF THE
                // PROCESSED LEVELS.
                // if (rf.bf.lambda != 0 && pivot_order->is_above_lambda(rf.bf.lambda, res_level)) {
                //     assert(down.first == 0);
                // }
                if (is_emptyset(irred_Ai))
                    break; // nothing left to be reduced.
            }
        }
        C_irreducibles1->setFull(ZtoNode(a_val), irred_Ai);
    }

    // cleanup
    MEDDLY::unpacked_node::Recycle(B);
    MEDDLY::unpacked_node::Recycle(A);
    // reduce and return result
    result.first  = res1F->createReducedNode(-1, C_irreducibles1);
    result.second = (compute_differences ? res2F->createReducedNode(-1, C_reduced2) : 0);
    saveResult(key, /*a, divisor,*/ result);

    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

reduce* 
reduce_opname::buildOperation(MEDDLY::forest *forestMDD,
                              const variable_order *pivot_order,
                              const bool compute_differences)
{
    if (0==forestMDD) return 0;

    if (forestMDD->isForRelations())
        throw MEDDLY::error(MEDDLY::error::TYPE_MISMATCH, __FILE__, __LINE__);

    if (forestMDD->getEdgeLabeling() == MEDDLY::edge_labeling::MULTI_TERMINAL) {
        if (forestMDD->isForRelations())
            throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED);

        return new reduce(this, forestMDD, pivot_order, compute_differences);
    }
    throw MEDDLY::error(MEDDLY::error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

/////////////////////////////////////////////////////////////////////////////////////////
















// /////////////////////////////////////////////////////////////////////////////////////////
// // Completion Procedure
// // Remove elements that are equal between [0..j-1] and smaller at var j
// /////////////////////////////////////////////////////////////////////////////////////////

// compl_proc::compl_proc(MEDDLY::binary_opname* opcode, 
// /**/                   MEDDLY::forest* arg1, MEDDLY::forest* arg2, 
// /**/                   MEDDLY::forest* res,
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

//     for (size_t i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
//         if (a_full && 0==A->down(i))
//             continue;
//         int a_val = NodeToZ(a_full ? i : A->index(i));

//         for (size_t j = 0; j < (b_full ? b_size : B->getSize()); j++) { // for each b
//             if (b_full && 0==B->down(j))
//                 continue;
//             int b_val = NodeToZ(b_full ? j : B->index(j));

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

//                 MEDDLY::node_handle res_ij = compute(A->down(i), B->down(j));
//                 unionNodes(C, res_ij, idx, resF, mddUnion);
//             }
//         }
//     }

//     // cleanup
//     MEDDLY::unpacked_node::Recycle(B);
//     MEDDLY::unpacked_node::Recycle(A);
//     // reduce and return result
//     result = resF->createReducedNode(-1, C);
//     saveResult(key, a, b, result);
//     return result;
// }

// /////////////////////////////////////////////////////////////////////////////////////////

// compl_proc* 
// compl_proc_opname::buildOperation(MEDDLY::forest* a1, 
//                                   MEDDLY::forest* a2, 
//                                   MEDDLY::forest* r,
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
// compl_proc_opname::buildOperation(MEDDLY::forest* a1, 
//                                           MEDDLY::forest* a2, 
//                                           MEDDLY::forest* r)
// {
//     throw std::exception(); // Unimplemented
// }

// /////////////////////////////////////////////////////////////////////////////////////////

// compl_proc_table::compl_proc_table(MEDDLY::forest* forest) {
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
sign_canon_mdd_opname::buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF)
{
    if (0==arF || 0==resF) return 0;
    
    return new sign_canon_mdd_op(this, arF, resF);
}

sign_canon_mdd_op::~sign_canon_mdd_op() { }

/////////////////////////////////////////////////////////////////////////////////////////

sign_canon_mdd_op::sign_canon_mdd_op(MEDDLY::unary_opname* opcode, MEDDLY::forest* _argF,
                                     MEDDLY::forest* _resF)
/**/ : base_NtoN(opcode, _argF, _resF)
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
    const int a_level = argF->getNodeLevel(a);
    const size_t a_size = get_node_size(A);

    const int res_size = max(a_size, (size_t)ZtoNode(-NodeToZ(a_size-1)) + 1);
    check_level_bound(resF, a_level, res_size); // resize level if needed
    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, a_level, res_size);

    const bool a_full = A->isFull();

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        if (a_val > 0) { // keep as-is, end recursive descent
            unionNodes(C, resF->linkNode(A->down(i)), 
                        ZtoNode(a_val), resF, mddUnion);
        }
        else if (a_val == 0) { 
            // following a 0-path, not found yet the first nnz
            // continue descending
            unionNodes(C, compute(A->down(i)), 
                        ZtoNode(a_val), resF, mddUnion);
        }
        else { // invert vectors from here down. Use VMULT
            assert(ZtoNode(-a_val) < resF->getDomain()->getVariableBound(a_level));
            unionNodes(C, VMULT->compute(A->down(i), -1), 
                        ZtoNode(-a_val), resF, mddUnion);
        }
    }

    // // Cleanup
    MEDDLY::unpacked_node::Recycle(A);
    // Save in cache
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, result);
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////














/////////////////////////////////////////////////////////////////////////////////////////
// Multiply vectors by a constant
/////////////////////////////////////////////////////////////////////////////////////////

vmult_opname::vmult_opname() 
/**/ : opname("VectorMultiply")
{ }

MEDDLY::operation* 
vmult_opname::buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF)
{
    if (0==arF || 0==resF) return 0;
    
    return new vmult_op(this, arF, resF);
}

vmult_op::~vmult_op() { }

/////////////////////////////////////////////////////////////////////////////////////////

vmult_op::vmult_op(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
                   MEDDLY::forest* _resF)
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
    MEDDLY::unpacked_node *A = MEDDLY::unpacked_node::New(argF);
    argF->unpackNode(A, a, MEDDLY::FULL_OR_SPARSE);
    const int a_level = argF->getNodeLevel(a);
    const size_t a_size = get_node_size(A);

    const int res_size = a_size * abs(multiplier) + (multiplier<0 ? 1 : 0);
    check_level_bound(resF, a_level, res_size);
    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, a_level, res_size);

    const bool a_full = A->isFull();

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        int r_val = ZtoNode(multiply_exact(multiplier, a_val));

        unionNodes(C, compute(A->down(i), multiplier), 
                   r_val, resF, mddUnion);
    }

    // // Cleanup
    MEDDLY::unpacked_node::Recycle(A);
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

vdivide_mdd_opname::vdivide_mdd_opname() 
/**/ : opname("VecCanonicalize")
{ }

MEDDLY::operation* 
vdivide_mdd_opname::buildOperation(MEDDLY::forest* arF, MEDDLY::forest* resF)
{
    if (0==arF || 0==resF) return 0;
    
    return new vdivide_mdd_op(this, arF, resF);
}

vdivide_mdd_op::~vdivide_mdd_op() { }

/////////////////////////////////////////////////////////////////////////////////////////

vdivide_mdd_op::vdivide_mdd_op(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
                               MEDDLY::forest* _resF)
/**/ : base_NItoNN(opcode, _argF, _resF, _resF)
{ }

/////////////////////////////////////////////////////////////////////////////////////////

// <first=undivisible, second=divided by d>
std::pair<MEDDLY::node_handle, MEDDLY::node_handle>
vdivide_mdd_op::compute(MEDDLY::node_handle a, const int divisor) 
{
    if (a==0 || a==-1) return std::make_pair(0, a);
    if (divisor==1)    return std::make_pair(res1F->linkNode(a), 0);

    std::pair<MEDDLY::node_handle, MEDDLY::node_handle> result;
    MEDDLY::ct_entry_key* CTsrch;
    MEDDLY::ct_entry_key* key = findResult(a, divisor, result);
    if (nullptr==key)
        return result;

    // Read the node and accumulate the LCM of the GCDs
    MEDDLY::unpacked_node* A = argF->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    const int a_level = argF->getNodeLevel(a);
    const size_t a_size = get_node_size(A);

    const int res_size = a_size;
    MEDDLY::unpacked_node* Cn = MEDDLY::unpacked_node::newFull(res1F, a_level, res_size); // not divisible by @divisor
    MEDDLY::unpacked_node* Cy = MEDDLY::unpacked_node::newFull(res2F, a_level, res_size); // divided by @divisor

    const bool a_full = A->isFull();

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        if (0 == (abs(a_val) % divisor)) { // path is divisible by @divisor
            std::pair<MEDDLY::node_handle, MEDDLY::node_handle> down = compute(A->down(i), divisor);
            unionNodes(Cn, down.first, ZtoNode(a_val), res1F, mddUnion);
            unionNodes(Cy, down.second, ZtoNode(a_val / divisor), res2F, mddUnion);
        }
        else { // path is not divisible by @divisor -> copy/ref A->down(i)
            unionNodes(Cn, res1F->linkNode(A->down(i)), ZtoNode(a_val), res1F, mddUnion);
        }
    }

    // // Cleanup
    MEDDLY::unpacked_node::Recycle(A);
    // Save in cache
    result.first  = res1F->createReducedNode(-1, Cn);
    result.second = res2F->createReducedNode(-1, Cy);
    saveResult(key, /*a, divisor,*/ result);

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
    
    return new divisors_finder_mdd_op(this, arg);
}

/////////////////////////////////////////////////////////////////////////////////////////

divisors_finder_mdd_op::divisors_finder_mdd_op(MEDDLY::unary_opname* oc, MEDDLY::forest* arg) 
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

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        lut_key lk_i = compute(A->down(i));
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
    MEDDLY::unpacked_node::Recycle(A);
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
/**/        MEDDLY::forest* arg1, MEDDLY::forest* arg2, MEDDLY::forest* res,
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

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        for (unsigned j = 0; j < (b_full ? b_size : B->getSize()); j++) { // for each b
            if (b_full && 0==B->down(j))
                continue;
            int b_val = NodeToZ(b_full ? j : B->index(j));

            MEDDLY::node_handle n = 0; // n = a if supp(a) subseteq_neq supp(b)
            if (lambda != 0 && p_table->pivot_order->is_above_lambda(lambda, res_level)) {
                n = compute(A->down(i), B->down(j)); // no check
            }
            else if (a_val == 0) { // a==0 -> b can only be 0
                if (b_val == 0) {
                    n = compute(A->down(i), B->down(j));
                }
            }
            else { // a_val != 0 -> b can be any value or 0
                if (is_pot_eq_supp && (b_val==0 || abs(b_val) < abs(a_val))) 
                    // b has smaller support than a, or b is smaller than a (for reduction)
                    n = p_table->get_op(lambda, false, subtract)->compute(A->down(i), B->down(j));
                else //if (multiply_exact(a_val, b_val) >= 0)
                    n = compute(A->down(i), B->down(j));
            }
            if (n != 0)
                unionNodes(C, n, ZtoNode(subtract ? a_val - b_val : a_val), 
                           resF, mddUnion);
        }
    }

    // cleanup
    MEDDLY::unpacked_node::Recycle(B);
    MEDDLY::unpacked_node::Recycle(A);
    // reduce and return result
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, b, result);
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

support_inclusion_op* 
support_inclusion_opname::buildOperation(MEDDLY::forest* a1, 
                                         MEDDLY::forest* a2, 
                                         MEDDLY::forest* r,
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
support_inclusion_opname::buildOperation(MEDDLY::forest* a1, 
        MEDDLY::forest* a2, MEDDLY::forest* r)
{
    throw std::exception(); // Unimplemented
}

/////////////////////////////////////////////////////////////////////////////////////////

support_inclusion_table::support_inclusion_table(MEDDLY::forest* forest,
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
    
    return new lv_selector_op(this, arF, resF, lambda);
}

lv_selector_op::~lv_selector_op() { }

/////////////////////////////////////////////////////////////////////////////////////////

lv_selector_op::lv_selector_op(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
                             MEDDLY::forest* _resF, const size_t lambda)
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

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        // cout << "a_val:"<<a_val<<" a_level:"<<a_level<<" value:"<<value<<" lambda:"<<lambda<<endl;
        if (a_level != lambda || a_val == value) // select this path that has @value at level @lambda
            unionNodes(C, compute(A->down(i), value), ZtoNode(a_val), resF, mddUnion);
    }

    // Cleanup
    MEDDLY::unpacked_node::Recycle(A);
    // Save in cache
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, value, result);
    return result;
 }

/////////////////////////////////////////////////////////////////////////////////////////

lv_selector_table::lv_selector_table(MEDDLY::forest* forest)
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

void domain_values_enumerator::visit(const MEDDLY::forest *argF,
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

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

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
            visit(argF, A->down(i));
        }
    }

    // Cleanup
    MEDDLY::unpacked_node::Recycle(A);
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
    
    return new smallest_degree_op(this, arg, p_table, degtype, lambda);
}

/////////////////////////////////////////////////////////////////////////////////////////

smallest_degree_op::smallest_degree_op(MEDDLY::unary_opname* oc, 
                                       MEDDLY::forest* arg,
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

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        int down_degree = compute(A->down(i));

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
    MEDDLY::unpacked_node::Recycle(A);
    // Save in cache
    res = smallest_degree;
    saveResult(key, a, res);
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////

smallest_degree_table::smallest_degree_table(MEDDLY::forest* forest,
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
    
    return new degree_finder_op(this, arg, p_table, degtype, lambda);
}

/////////////////////////////////////////////////////////////////////////////////////////

degree_finder_op::degree_finder_op(MEDDLY::unary_opname* oc, 
                                       MEDDLY::forest* arg,
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

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));

        lut_key down_lk = compute(A->down(i));
        const std::vector<int>& down_degrees = p_table->look_up(down_lk);

        int edge_incr = 0;
        if (p_table->pivot_order->is_below_lambda(lambda, a_level)) { 
            // level contributes to the degree
            edge_incr = get_degree_of(a_val, degtype);
            assert(edge_incr >= 0);
        }

        for (int d : down_degrees)
            degrees.push_back(add_exact(d, edge_incr));
    }

    // Sort and make the list of degrees unique
    std::sort(degrees.begin(), degrees.end(), [](int i, int j){ return i>j;});
    degrees.erase(std::unique(degrees.begin(), degrees.end()), degrees.end());

    // Cleanup
    MEDDLY::unpacked_node::Recycle(A);
    // Save in cache
    res = p_table->lut.insert(std::move(degrees));
    saveResult(key, a, res);
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////

degree_finder_table::degree_finder_table(MEDDLY::forest* forest,
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
    
    return new degree_selector_op(this, arF, resF, p_table, degtype, lambda);
}

degree_selector_op::~degree_selector_op() { }

/////////////////////////////////////////////////////////////////////////////////////////

degree_selector_op::degree_selector_op(MEDDLY::opname* opcode, MEDDLY::forest* _argF,
                                       MEDDLY::forest* _resF,
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

    MEDDLY::unpacked_node* A = argF->newUnpacked(a, MEDDLY::SPARSE_ONLY);
    const int a_level = argF->getNodeLevel(a);
    const size_t a_size = get_node_size(A);

    const int res_size = a_size;
    check_level_bound(resF, a_level, res_size);
    MEDDLY::unpacked_node* C = MEDDLY::unpacked_node::newFull(resF, a_level, res_size);

    const bool a_full = A->isFull();

    for (unsigned i = 0; i < (a_full ? a_size : A->getSize()); i++) { // for each a
        if (a_full && 0==A->down(i))
            continue;
        int a_val = NodeToZ(a_full ? i : A->index(i));
        
        int r_val = a_val;
        int down_degree = degree;
        if (p_table->pivot_order->is_below_lambda(lambda, a_level)) { 
            // level contributes to the degree
            down_degree -= get_degree_of(a_val, degtype);
        }

        unionNodes(C, compute(A->down(i), down_degree), 
                   ZtoNode(r_val), resF, mddUnion);
    }

    // // Cleanup
    MEDDLY::unpacked_node::Recycle(A);
    // Save in cache
    result = resF->createReducedNode(-1, C);
    saveResult(key, a, degree, result);
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////

degree_selector_table::degree_selector_table(MEDDLY::forest* forest,
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
    // Declare out custom operators
    S_VECTORS_OPNAME = new s_vectors_opname();
    S_VECTORS = new s_vectors(S_VECTORS_OPNAME, forestMDD, forestMDD, forestMDD, pivot_order);

    // LESSEQ_SQ_OPNAME = new lesseq_sq_opname();
    // LESSEQ_SQ_OPS = new lesseq_sq_table(forestMDD, pivot_order);

    LEQ_NEQ_SQ_OPNAME = new leq_neq_sq_opname();
    LEQ_NEQ_SQ_COMPARE = LEQ_NEQ_SQ_OPNAME->buildOperation(forestMDD, pivot_order, false);
    LEQ_NEQ_SQ_SUBTRACT = LEQ_NEQ_SQ_OPNAME->buildOperation(forestMDD, pivot_order, true);

    REDUCE_OPNAME = new reduce_opname();
    REDUCE = REDUCE_OPNAME->buildOperation(forestMDD, pivot_order, true);
    GET_IRREDUCIBLES = REDUCE_OPNAME->buildOperation(forestMDD, pivot_order, false);

    // REDUCE3_OPNAME = new reduce3_opname();
    // REDUCE3 = REDUCE3_OPNAME->buildOperation(forestMDD, pivot_order);

    // COMPL_PROC_OPNAME = new compl_proc_opname();
    // COMPL_PROC_OPS = new compl_proc_table(forestMDD);

    SIGN_CANON_OPNAME = new sign_canon_mdd_opname();
    SIGN_CANON = SIGN_CANON_OPNAME->buildOperation(forestMDD, forestMDD);

    VMULT_OPNAME = new vmult_opname();
    VMULT = (vmult_op*)VMULT_OPNAME->buildOperation(forestMDD, forestMDD);

    VDIVIDE_OPNAME = new vdivide_mdd_opname();
    VDIVIDE = (vdivide_mdd_op*)VDIVIDE_OPNAME->buildOperation(forestMDD, forestMDD);

    DIV_FINDER_MDD_OPNAME = new divisors_finder_mdd_opname();
    DIV_FINDER_MDD = (divisors_finder_mdd_op*)DIV_FINDER_MDD_OPNAME->buildOperation(forestMDD);

    SUPPORT_INCL_OPNAME = new support_inclusion_opname();
    SUPPORT_INCL_TABLE = new support_inclusion_table(forestMDD, pivot_order);

    LV_SELECTOR_OPNAME = new lv_selector_opname();
    LV_SELECTOR_OPS = new lv_selector_table(forestMDD);

    SMALLEST_DEGREE_OPNAME = new smallest_degree_opname();
    SMALLEST_DEGREE_TABLE = new smallest_degree_table(forestMDD, pivot_order);

    DEGREE_FINDER_OPNAME = new degree_finder_opname();
    DEGREE_FINDER_TABLE = new degree_finder_table(forestMDD, pivot_order);

    DEGREE_SELECTOR_OPNAME = new degree_selector_opname();
    DEGREE_SELECTOR_TABLE = new degree_selector_table(forestMDD, pivot_order);
}

/////////////////////////////////////////////////////////////////////////////////////////

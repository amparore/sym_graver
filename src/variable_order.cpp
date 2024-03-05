#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cassert>
#include <set>
#include <array>
#include <stack>
#include <cstdlib>

#include "variable_order.h"
#include "matrix.h"

#include "reverse_heap.h"

#ifdef HAS_BOOST_CPP
// Boost graph 
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/sloan_ordering.hpp>
#endif // HAS_BOOST_CPP

using namespace std;


/////////////////////////////////////////////////////////////////////////////////////////

variable_order::variable_order(size_t m) {
    var_to_level.resize(m);
    level_to_var.resize(m);
    std::iota(var_to_level.begin(), var_to_level.end(), 0);
    std::iota(level_to_var.begin(), level_to_var.end(), 0);
}

/////////////////////////////////////////////////////////////////////////////////////////

void variable_order::check_order() const {
    std::vector<bool> has_var(var_to_level.size());
    for (size_t i=0; i<var_to_level.size(); i++) {
        if (has_var[lvl2var(i)])
            throw runtime_error("Incorrect variable order: repetition.");
        has_var[lvl2var(i)] = true;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

void variable_order::print() const {
    cout << "  lvl2var: ";
    for (size_t lvl=0; lvl<level_to_var.size(); lvl++)
        cout << (lvl2var(lvl) + 1) << " ";
    cout << endl;

    cout << "  var2lvl: ";
    for (size_t var=0; var<level_to_var.size(); var++)
        cout << (var2lvl(var) + 1) << " ";
    cout << endl << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int>>
reorder_matrix(const std::vector<std::vector<int>>& A,
               const variable_order &vorder)
{
    std::vector<std::vector<int>> B = A;
    const size_t n = A.size();
    if (n > 0) {
        const size_t m = A[0].size();
        for (size_t j1=0; j1<m; j1++) {
            size_t j2 = vorder.var2lvl(j1);
            for (size_t i=0; i<n; i++) // copy column
                B[i][j2] = A[i][j1];
        }
    }
    return B;
}

/////////////////////////////////////////////////////////////////////////////////////////

// void canonicalize_by_order(std::vector<int>& row,
//                            const variable_order &vorder)
// {
//     size_t j=0;
//     // search first non-zero
//     while (j<row.size()) {
//         size_t k = vorder.lvl2var(j);
//         cout << "i="<<0<<" j="<<j<<" k="<<k<<"  row[k]="<<row[k]<<endl;
//         if (row[k] > 0)
//             return; // positive: row is canonical
//         else if (row[k] < 0)
//             break; // negative
//         j++;
//     }
//     // invert sign of @row
//     while (j<row.size()) {
//         size_t k = vorder.lvl2var(j);
//         row[k] = -row[k];
//         j++;
//     }
// }

// void canonicalize_by_order(std::vector<std::vector<int>>& mat,
//                            const variable_order &vorder)
// {
//     for (size_t i=0; i<mat.size(); i++) {
//         canonicalize_by_order(mat[i], vorder);
//     }
// }

// /////////////////////////////////////////////////////////////////////////////////////////











/////////////////////////////////////////////////////////////////////////////////////////

class mat_graph {
public:
    // size of the graph
    const size_t m;
    // incidence
    std::vector<std::vector<size_t>> inc;
    
    mat_graph(const std::vector<std::vector<int>>& A);

    void save_dot(const char* fname) const;
    // out-degree of vertex j
    inline size_t out_degree_of(size_t j) const { return inc[j].size(); }
};

ostream& operator<<(ostream&, const mat_graph& mg);

/////////////////////////////////////////////////////////////////////////////////////////

mat_graph::mat_graph(const std::vector<std::vector<int>>& A) : m(A[0].size()), inc(m)
{
    std::set<pair<int, int>> inserted;

    // Connect variables that share the same constraint
    std::vector<size_t> vars_in_constr;
    for (size_t i=0; i<A.size(); i++) { // for each constraint
        vars_in_constr.clear();
        for (size_t j=0; j<A[i].size(); j++) {
            if (A[i][j] != 0)
                vars_in_constr.push_back(j);
        }

        // add edges
        for (size_t j1=0; j1<vars_in_constr.size(); j1++) {
            for (size_t j2=j1+1; j2<vars_in_constr.size(); j2++) {
                const size_t v1 = vars_in_constr[j1];
                const size_t v2 = vars_in_constr[j2];
                auto key = std::make_pair(v1 < v2 ? v1 : v2,
                                          v1 < v2 ? v2 : v1);
                if (inserted.count(key) == 0) {
                    inc[v1].push_back(v2);
                    inc[v2].push_back(v1);
                    inserted.insert(key);
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<<(ostream& os, const mat_graph& mg) {
    os << "Graph G("<<mg.m<<");\n";
    for (size_t i=0; i<mg.m; i++) {
        os << "  From "<<i<<":";
        for (size_t j : mg.inc[i])
            os << " " << j;
        os << endl;
    }
    os << endl;
    return os;
}

/////////////////////////////////////////////////////////////////////////////////////////

void mat_graph::save_dot(const char* fname) const {
	ofstream og(fname);
	og << "graph { " << endl;
    for (size_t i=0; i<m; i++) {
        for (size_t j : inc[i]) {
            if (i < j) // avoid duplicating bidirectional edges.
                og << "v" << i << " -- v" << j << ";\n";
        }
    }
	og << "}\n" << endl;
	og.close();
}

/////////////////////////////////////////////////////////////////////////////////////////

class incremental_reorderer_base {
protected:
    const mat_graph& graph;  // incidence graph
    std::vector<bool> S;    // already selected variables
    // Heap of the variables, in weight order
    reverse_heap<double> RH;

    virtual double get_variable_weight(size_t i) = 0;
    virtual void select_variable(size_t i, size_t seq_num) = 0;
    virtual void update_end() = 0;
public:
    incremental_reorderer_base(const mat_graph& g) : graph(g), RH(g.m) { }

    void compute(variable_order &out_order) {
        S.resize(graph.m, false);
        // std::vector<int> update_list;
        // update_list.reserve(graph.m);

        // Initialize all weights and insert them into the heap
        for (int var=0; var<graph.m; var++)
            RH.push_heap(var, get_variable_weight(var));

        // Generate the variable order
        for (int pos=0; pos<graph.m; pos++) {
            double top_w = RH.top_weight();
            int var = RH.pop_heap(); // variable with the highest weight
            // cout << ""<<var<<"  with weight "<<top_w<<endl;
            S[var] = true;
            out_order.bind_var2lvl(var, pos);
            select_variable(var, pos);

            // Update the weights
            for (auto upd_var : graph.inc[var]) {
                if (!S[upd_var])
                    RH.update_weight(upd_var, get_variable_weight(upd_var));
            }
            update_end();
            // update_list.resize(0);
        }
    }
};

enum vertex_status {
    POSTACTIVE, ACTIVE, PREACTIVE, INACTIVE
};

class incremental_reorderer_mininc : public incremental_reorderer_base {
public:
    incremental_reorderer_mininc(const mat_graph& g, double _w1, double _w2) 
    /**/ : incremental_reorderer_base(g), w1(_w1), w2(_w2) {
        status.resize(graph.m, INACTIVE);
        var_boost.resize(graph.m, 0.0);
    }

protected:
    double w1, w2;
    std::vector<vertex_status> status;
    std::vector<int> var_boost;

    virtual void select_variable(size_t i, size_t seq_num) override {
        assert(status[i]==INACTIVE || status[i]==ACTIVE);
        // set i as postactive, and its neighbours as active
        status[i] = POSTACTIVE;
        for (size_t adj : graph.inc[i]) {
            if (status[adj] == INACTIVE) {//(!S[adj])
                status[adj] = ACTIVE;
                var_boost[adj] += seq_num; // Sloan-like boost
            }
        }
    }

    virtual double get_variable_weight(size_t i) override {
        // count the number of non-selected adjacent variables
        size_t actives = 0, preactives = 0, selected = 0;
        for (size_t adj : graph.inc[i]) {
            if (!S[adj]) {
                if (status[adj] == ACTIVE)
                    ++actives;
                else 
                    ++preactives;
            }
            else 
                ++selected;
        }

        const size_t degree = graph.inc[i].size();

        double w = w1*double(preactives) + actives + w2*selected + var_boost[i];
        w /= degree + 0.1;

        if (status[i]==INACTIVE)
            w -= 100; // do not favor inactive nodes

        // cout << "get_variable_weight("<<i<<") = "<<w
        //      << "  preactives="<<preactives<<" actives="<<actives<<" degree="<<degree<<endl;
        return w;
    }

    void update_end() override {
        // for (size_t i=0; i<graph.m; i++) {
        //     double w = RH.get_weight(i);
        //     if (0==w)
        //         cout << "   *";
        //     else
        //         cout << setw(4) << setprecision(1) << w;
        //     const char *s[4] = {"-", "A", "P", "I"};
        //     cout << s[status[i]] << " ";
        // }
        // cout << endl;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////

void fast_varorder(const std::vector<std::vector<int>>& A,
                   variable_order &out_order)
{
    mat_graph mg(A);
    // cout << mg << endl;
    // mg.save_dot("varorder.dot");

    // std::array<std::array<double, 2>, 4> weigths = {{ 
    //     {-2.0, 4.0}, {-3.0, 3.0}, {-6.0, 1.0}, {-2.0, 8.0} 
    // }};

    // size_t best_score = size_t(-1);
    // for (const auto& w : weigths) {
    //     variable_order order(mg.m);
    //     incremental_reorderer_mininc reordered(mg, w[0], w[1]);
    //     reordered.compute(order);

    //     auto rA = A;
    //     reorder_matrix(rA, order);
    //     size_t score = irank(rA);
    //     cout << "irank = " << score << endl;
    //     if (best_score == size_t(-1) || best_score > score) {
    //         best_score = score;
    //         out_order = order;
    //     }
    // }

    incremental_reorderer_mininc reordered(mg, -4.0, 2.0);
    reordered.compute(out_order);
}

/////////////////////////////////////////////////////////////////////////////////////////








/////////////////////////////////////////////////////////////////////////////////////////
// Sloan variable reordering algorithm
// Described in:
//  Sloan, S. W. "A FORTRAN program for profile and wavefront reduction." 
//  International Journal for Numerical Methods in Engineering 28.11 (1989): 2651-2679.
/////////////////////////////////////////////////////////////////////////////////////////

class root_level_structure {
    const mat_graph& mg; // visited graph
    std::vector<size_t> enqueued; // visited state
    std::stack<std::pair<size_t, size_t>> q; // <vertex, depth> visit stack

public:
    inline root_level_structure(const mat_graph& mg) 
    /**/ : mg(mg), enqueued(mg.m, size_t(-1)) {}

    // visit the graph starting from v0, and compute the maximum width and the total depth
    // Optionally, the distance from v0 can be stored in p_distances[]
    void dfs_width_depth(const size_t v0, size_t &max_width, size_t &max_depth,
                         std::vector<size_t> *p_distances) 
    {
        max_width = max_depth = 0;
        assert(q.empty());
        q.push(make_pair(v0, 0));
        enqueued[v0] = v0;
        size_t d_prev = 0, curr_width = 0;
        // do the DFS visit
        while (!q.empty()) {
            size_t v, d;
            std::tie(v, d) = q.top();
            q.pop();
            if (p_distances)
                (*p_distances)[v] = d;
            max_depth = std::max(max_depth, d);
            ++curr_width;
            if (d_prev != d) { // visiting new depth layer in DFS order
                max_width = std::max(max_width, curr_width);
                curr_width = 0;
                d_prev = d;
            }
            for (size_t adj : mg.inc[v]) {
                if (v0 != enqueued[adj]) {
                    q.push(make_pair(adj, d+1));
                    enqueued[adj] = v0;
                }
            }
        }
        max_width = std::max(max_width, curr_width);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////

void sloan_varorder(const mat_graph& mg,
                    variable_order &out_order) 
{
    double W1=1, W2=2;
    // get the root-level-structure scores for each node in the graph
    root_level_structure rls(mg);
    struct wd_t { size_t width, depth; };
    std::vector<wd_t> rlss(mg.m); // root level structure scores
    for (size_t i=0; i<mg.m; i++) {
        rls.dfs_width_depth(i, rlss[i].width, rlss[i].depth, nullptr);
        // cout << "  vertex:"<<i<<"   width:"<<rlss[i].width<<" depth:"<<rlss[i].depth<<endl;
    }
    
    std::vector<bool> visited(mg.m, false);
    std::vector<size_t> distances(mg.m, size_t(-1));
    std::vector<vertex_status> status(mg.m, INACTIVE);
    reverse_heap<double> RH(mg.m);
    size_t order_position = 0;

    // Get the vertex with the smallest width/highest depth among the ones not yet visited
    while (order_position < mg.m) 
    {
        size_t start = size_t(-1);
        for (size_t i=0; i<mg.m; i++) {
            if (visited[i])
                continue;

            if (start==size_t(-1) || 
                ((rlss[i].depth > rlss[start].depth) && 
                 (rlss[i].width < rlss[start].width)))
            {
                start = i;
            }
        }
        assert(start != size_t(-1));

        // Now find the matching end vertex that keeps the smallest tree width
        size_t sd, sw, end=size_t(-1), ew, ed;
        rls.dfs_width_depth(start, sw, sd, &distances);
        // cout << "start "<<start<<" has sd="<<sd<<endl;
        for (size_t i=0; i<mg.m; i++) {
            if (visited[i])
                continue;
            // cout << "d["<<i<<"]="<<distances[i]<<" ";
            if (distances[i] == sd) { // i is at maximum distance from start
                if (end==size_t(-1) ||
                    (rlss[i].width < rlss[end].width))
                {
                    end = i;
                }
            }
        }

        // cout << "start="<<start<<" end="<<end<<endl;

        // Now visit from the end, and store the distances
        rls.dfs_width_depth(end, ew, ed, &distances);
        
        // priority-ordered visit from start. 
        // Status of vertices guide the steps (following Sloan algorithm).
        //  INACTIVE    Not yet inserted in the priority queue
        //  PREACTIVE   Inserted in the priority queue, not yet visited its neighbours
        //  ACTIVE      Vertex and its neighbours are inserted in the priority queue
        //  POSTACTIVE  Visit completed.
        RH.push_heap(start, 1.0);
        status[start] = PREACTIVE;

        while (!RH.empty()) {
            size_t i = RH.pop_heap();
            assert(status[i] != POSTACTIVE);

            // Update adjacent (depth=1 visit)
            if (status[i] == PREACTIVE) { // activate j
                for (size_t j : mg.inc[i]) { 
                    if (status[j] == INACTIVE) { // insert j
                        status[j] = PREACTIVE;
                        RH.push_heap(j, W1*distances[j] - W2*(mg.out_degree_of(j)+1));
                    }
                    else if (status[j] != POSTACTIVE) { // boost j
                        RH.update_weight(j, RH.get_weight(j) + W2);
                    }
                }
            }

            status[i] = POSTACTIVE;
            visited[i] = true;
            out_order.bind_var2lvl(i, order_position++);
            // cout << "visiting "<<i<<endl;

            // Pre-activate the nodes at distance 2
            for (size_t j : mg.inc[i]) { 
                if (status[j] == PREACTIVE) { // activate j
                    status[j] = ACTIVE;
                    RH.update_weight(j, RH.get_weight(j) + W2);

                    // visit neighbours of j and pre-activate them
                    for (size_t k : mg.inc[j]) {
                        if (status[k] == INACTIVE) { // insert k
                            status[k] = PREACTIVE;
                            RH.push_heap(k, W1*distances[k] - W2*(mg.out_degree_of(k)+1));
                        }
                        else if (status[k] != POSTACTIVE) { // boost k
                            RH.update_weight(k, RH.get_weight(k) + W2);
                        }
                    }
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

void sloan_varorder(const std::vector<std::vector<int>>& A,
                     variable_order &out_order)
{
    mat_graph mg(A);
    sloan_varorder(mg, out_order);
}

/////////////////////////////////////////////////////////////////////////////////////////






















#ifdef HAS_BOOST_CPP

/////////////////////////////////////////////////////////////////////////////////////////
// Sloan variable reordering
/////////////////////////////////////////////////////////////////////////////////////////

template<typename Graph>
void init_graph_edges_from_matrix(Graph &graph, 
                                  const std::vector<std::vector<int>>& A)
{
    std::set<pair<int, int>> inserted;

    // Connect variables that share the same constraint
    std::vector<size_t> vars_in_constr;
    for (size_t i=0; i<A.size(); i++) { // for each constraint
        vars_in_constr.clear();
        for (size_t j=0; j<A[i].size(); j++) {
            if (A[i][j] != 0)
                vars_in_constr.push_back(j);
        }

        // add edges
        for (size_t j1=0; j1<vars_in_constr.size(); j1++) {
            for (size_t j2=j1+1; j2<vars_in_constr.size(); j2++) {
                const size_t v1 = vars_in_constr[j1];
                const size_t v2 = vars_in_constr[j2];
                auto key = std::make_pair(v1 < v2 ? v1 : v2,
                                          v1 < v2 ? v2 : v1);
                if (inserted.count(key) == 0) {
                    boost::add_edge(v1, v2,  graph);
                    boost::add_edge(v2, v1,  graph);
                    inserted.insert(key);
                }
            }
        }
    }

	// ofstream og("graph.dot");
	// og << "graph { " << endl;
	// typedef typename boost::graph_traits<Graph>::edge_iterator  edge_iterator;
	// edge_iterator e, end;
	// for (boost::tie(e, end) = boost::edges(graph); e != end; ++e) {
	// 	og << "v" << boost::source(*e, graph) << " -- v" << boost::target(*e, graph) << ";\n";
	// }
	// og << "}\n" << endl;
	// og.close();

    // cout << "Graph G("<<boost::num_vertices(graph)<<");\n";
    // typedef typename boost::graph_traits<Graph>::edge_iterator  edge_iterator;
    // edge_iterator e, end;
    // for (boost::tie(e, end) = boost::edges(graph); e != end; ++e) {
    //     cout << "boost::add_edge(" << boost::source(*e, graph) << ", " 
    //          << boost::target(*e, graph) << ", graph);\n";
    // }
}

/////////////////////////////////////////////////////////////////////////////////////////

// Convert a boost::graph permutation into a variable ordering for Meddly
template<typename VertexIndexMap, typename VectorPerm>
void fill_out_ordering(variable_order &out_order, 
                       const size_t num_vars,
                       const VectorPerm& permutation,
                       const VertexIndexMap& vertex_index_map) 
{
	typedef typename VectorPerm::size_type  size_type;
    // Print the ordering
    // cout << permutation.size() << " VALUES: ";
    // for (size_type c = 0; c != permutation.size(); ++c) {
    //     size_type index = vertex_index_map[ permutation[c] ];
    //     if (index >= graph.m) // the boost::graph vertex is not a place
    //         cout << "<> ";
    //     else
    //         cout << index << " ";
    //         // cout << tabp[ index ].place_name << " ";
    // }
    // cout << endl << endl;
    // for (size_type var=0; var<graph.m; var++) {
    // 	if (permutation.end() != std::find(permutation.begin(), permutation.end(), var))
    // 		continue;
    // 	cout << "missing variable " << var << " in ordering..." << endl;
    // }

    // Copy the ordering into out_order[]
    size_t level = 0;
    for (size_type c = 0; c != permutation.size(); ++c) {
        size_type index = vertex_index_map[ permutation[c] ];
        if (index >= num_vars) // the boost::graph vertex is not a place
            continue;
        if (index < 0)
            throw runtime_error("Negative index!");
        out_order.bind_var2lvl(index, level++);
        // out_order[ index ] = level++;
    }
    if (level != num_vars)
    	throw runtime_error("Ordering is incomplete.");
    out_order.check_order(); // check repetitions
}

/////////////////////////////////////////////////////////////////////////////////////////

void boost_sloan_varorder(const std::vector<std::vector<int>>& A,
                          variable_order &out_order) 
{
    const size_t m = A[0].size();
    // The graph type for the sloan method
    typedef boost::adjacency_list<boost::setS, 
    /**/                          boost::vecS, 
    /**/                          boost::undirectedS, 
    /**/                          boost::property<boost::vertex_color_t, 
    /**/                                          boost::default_color_type,
    /**/                          boost::property<boost::vertex_degree_t,
    /**/                                          int,
    /**/                          boost::property<boost::vertex_priority_t,
    /**/                                          double>>>> Graph;

    typedef boost::graph_traits<Graph>::vertex_descriptor  Vertex;
    typedef boost::graph_traits<Graph>::vertices_size_type size_type;

    Graph G(m);
    boost::property_map<Graph, boost::vertex_index_t>::type index_map = boost::get(boost::vertex_index, G);
    boost::property_map<Graph, boost::vertex_degree_t>::type degree_map = boost::get(boost::vertex_degree, G);
    init_graph_edges_from_matrix(G, A);

    // Setup vertex degrees
    boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    for (boost::tie(ui, ui_end) = boost::vertices(G); ui != ui_end; ++ui)
        degree_map[*ui] = boost::degree(*ui, G);

    // Support vector needed by sloan_ordering  
    std::vector<Vertex> sloan_order(boost::num_vertices(G));
    std::vector<size_type> perm(boost::num_vertices(G));
    std::vector<Vertex>::iterator inserted_vars;
    const int W1 = 1, W2 = 2;

    inserted_vars = boost::sloan_ordering(G, sloan_order.begin(), 
				                          boost::get(boost::vertex_color, G),
				                          boost::make_degree_map(G), 
				                          boost::get(boost::vertex_priority, G),
                                          W1, W2);

    if (inserted_vars != sloan_order.end()) {
    	// We missed some vertices. This may happen if those vertices are isolated.
    	// Complete the ordering by finding those vertices
    	typedef boost::property_map<Graph, boost::vertex_color_t>::type ColorMap;
	    typedef boost::property_traits<ColorMap>::value_type ColorValue;
	    typedef boost::color_traits<ColorValue> Color;
	    ColorMap color_map = boost::get(boost::vertex_color, G);
	    for (boost::tie(ui, ui_end) = boost::vertices(G); ui != ui_end; ++ui) {
    		if (color_map[*ui] == Color::white()) {
    			*inserted_vars++ = *ui;
    			color_map[*ui] = Color::green();
    		}
    	}
    }

    fill_out_ordering(out_order, m, sloan_order, index_map);
}

/////////////////////////////////////////////////////////////////////////////////////////

#endif // HAS_BOOST_CPP















/////////////////////////////////////////////////////////////////////////////////////////

// void pivot_order_from_matrix(variable_order& pivots,
//                              const std::vector<std::vector<int>>& A,
//                              bool optimize_graver) 
// {
//     const size_t n = A.size(), m = A.front().size();
//     // Count all positive/negative values on A columns
//     std::vector<std::pair<long, size_t>> weights(m);
//     for (size_t j=0; j<m; j++) {
//         int w;
//         if (optimize_graver) { // take interval / gcd
//             int max_pos = 0, max_neg = 0, g = 0;
//             for (size_t i=0; i<n; i++) {
//                 if (A[i][j] > 0) {
//                     max_pos = max(max_pos, A[i][j]);
//                     g = gcd(g, A[i][j]);
//                 }
//                 else if (A[i][j] < 0) {
//                     max_neg = max(max_neg, -A[i][j]);
//                     g = gcd(g, -A[i][j]);
//                 }
//             }
//             w = max_pos + max_neg;
//             if (g > 1)
//                 w /= g;
//         }
//         else { // take (sum(positives) * sum(negatives)) + sum
//             int pos=0, neg=0;
//             for (size_t i=0; i<n; i++) {
//                 if (A[i][j] > 0) 
//                     pos += A[i][j];
//                 else if (A[i][j] < 0) 
//                     neg -= A[i][j];
//             }
//             w = (pos * neg) + pos + neg;
//         }
//         weights[j] = make_pair(w, j);
//     }
//     std::sort(weights.begin(), weights.end());
//     for (size_t j=0; j<m; j++) {
//         // cout << "j:"<<j<<" w:"<<weights[j].first<<" var:"<<weights[j].second<<endl;
//         pivots.bind_var2lvl(j, weights[j].second); 
//         // pivots.bind_var2lvl(weights[j].second, j);
//     }
// }

/////////////////////////////////////////////////////////////////////////////////////////

void pivot_order_from_matrix_iter(variable_order& pivots,
                                  const std::vector<std::vector<int>>& A,
                                  const bool optimize_graver,
                                  const size_t num_iters,
                                  const std::vector<size_t>& fixed_vars) 
{
    const size_t n = A.size();
    if (n == 0)
        return;

    const size_t m = A.front().size();
    std::vector<bool> blocked_vars(m, false);
    for (size_t var : fixed_vars)
        blocked_vars[var] = true;

    typedef double weight_t;
    std::vector<weight_t> row_weights(n, weight_t(1));
    std::vector<weight_t> col_weights(m, weight_t(1));

    for (size_t step=0; step<num_iters; step++) {
        // Accumulate the row weights as the sum of all the involved column weights
        if (step > 0) {
            for (size_t i=0; i<n; i++) {
                weight_t rw = 0.0;
                for (size_t j=0; j<m; j++) {
                    if (!blocked_vars[j]) {
                        if (A[i][j] != 0)
                            rw += col_weights[j]; // column/var j is affected by row i
                    }
                }
                row_weights[i] = max(weight_t(1), rw);
            }
        }

        // recompute column/variable weights making an estimate of the total operation count
        for (size_t j=0; j<m; j++) {
            if (!blocked_vars[j]) {
                weight_t w;
                if (optimize_graver) { // take interval / gcd
                    weight_t max_pos = 0, max_neg = 0;
                    int g = 0;
                    for (size_t i=0; i<n; i++) {
                        if (A[i][j] > 0) {
                            max_pos = max(max_pos, A[i][j] * row_weights[i]);
                            g = gcd(g, A[i][j]);
                        }
                        else if (A[i][j] < 0) {
                            max_neg = max(max_neg, -A[i][j] * row_weights[i]);
                            g = gcd(g, -A[i][j]);
                        }
                    }
                    w = max_pos + max_neg;
                    if (g > 1)
                        w /= g;
                }
                else { // take (sum(positives) * sum(negatives)) + sum
                    weight_t pos=0, neg=0;
                    for (size_t i=0; i<n; i++) {
                        if (A[i][j] > 0) 
                            pos += A[i][j] * row_weights[i];
                        else if (A[i][j] < 0) 
                            neg -= A[i][j] * row_weights[i];
                    }
                    w = (pos * neg) + pos + neg;
                }
                col_weights[j] = w;
            }
        }

        // cout << "Step: "<<step << endl;
        // for (size_t i=0; i<n; i++) cout << row_weights[i] << " ";
        // cout << endl;
        // for (size_t j=0; j<m; j++) cout << (j+1)<<":"<< col_weights[j] << " ";
        // cout << endl;
    }

    // sort by ascending weight, then setup the pivot order
    std::vector<std::pair<weight_t, size_t>> sorted_weights(m);
    for (size_t j=0; j<m; j++)
        sorted_weights[j] = make_pair(col_weights[j], j);
    std::sort(sorted_weights.begin(), sorted_weights.end());
    // first put all fixed vars
    size_t pos = 0;
    for (size_t var : fixed_vars)
        pivots.bind_var2lvl(pos++, var); 
    // then put all the other vars, in weighted order
    for (size_t j=0; j<m; j++) {
        size_t var = sorted_weights[j].second;
        // cout << "j:"<<j<<" w:"<<sorted_weights[j].first<<" var:"<<sorted_weights[j].second<<endl;
        if (!blocked_vars[var])
            pivots.bind_var2lvl(pos++, var); 
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

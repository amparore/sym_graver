#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <algorithm>
#include <cassert>

#include "dd_operations.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////

// Draw a dd_edge in Graphviz format
struct dot_of_DD {
    ofstream dot;
    std::set<MEDDLY::node_handle> visited;
    MEDDLY::expert_forest *forest;
    std::vector<std::vector<std::string>> nodes_per_lvl; // unprimed and primed
    // const std::vector<bool> *pSingletonLevel;
    // const RSRG *rs;
    bool write_level_labels = false;
    bool write_terminals = false;

    void visit(const MEDDLY::node_handle node) {
        if (visited.count(node) > 0)
            return;
        if (node == forest->handleForValue(false) || 
            (node == forest->handleForValue(true) && !write_terminals))
            return;
        visited.insert(node);

        if (node == forest->handleForValue(true)) {
            dot << "  nT [label=\"T\"];\n";
            return;
        }

        const int node_level = forest->getNodeLevel(node);

        // MEDDLY::unpacked_node *rnode = MEDDLY::unpacked_node::newFromNode(forest, node, MEDDLY::unpacked_node::storage_style::FULL_NODE);
        MEDDLY::unpacked_node *rnode = MEDDLY::unpacked_node::New();
        forest->unpackNode(rnode, node, MEDDLY::FULL_ONLY);

        assert(rnode->isFull());
        int end = rnode->getSize() - 1; // find last non-false edge
        while (rnode->d(end) == forest->handleForValue(false))
            end--;

        int series = node_level<0 ? 1 : 0;
        // Specify the rank of the node
        nodes_per_lvl[series][abs(node_level)] += " n";
        nodes_per_lvl[series][abs(node_level)] += std::to_string(node);

        // draw the node
        ostringstream edges;
        dot << "  n"<<node<<" [label=\"";

        // determine sorting order
        std::vector<std::pair<int, int>> order;
        for (int i = 0; i <= end; i++) {
            if (rnode->d(i) == forest->handleForValue(false))
                continue;

            order.push_back(make_pair(NodeToZ(i), i));
        }
        std::sort(order.begin(), order.end());

        // write node entries (sorted by value)
        int cnt = 0;
        for (const auto& visit_pair : order) {
            int i = visit_pair.second;
            if (rnode->d(i) == forest->handleForValue(false))
                continue;
            
            dot << (cnt++==0 ? "" : "|") << "<i"<<i<<">" << NodeToZ(i);
            if (rnode->d(i) == forest->handleForValue(true)) {
                if (write_terminals)
                    edges << "  n"<<node<<":i"<<i<<" -> nT;\n";
            }
            else {
                edges << "  n"<<node<<":i"<<i<<" -> n"<<rnode->d(i)<<":n;\n";
            }
        }
        dot << "\"];\n";
        dot << edges.str();

        // Visit recursively
        for (int i = 0; i <= end; i++)
            visit(rnode->d(i));

        MEDDLY::unpacked_node::recycle(rnode);
    }

    void start_visit(const MEDDLY::dd_edge& e) {
        forest = static_cast<MEDDLY::expert_forest *>(e.getForest());
        bool isForRel = forest->isForRelations();
        int num_series = isForRel ? 2 : 1; // number of sets in nodes_per_lvl[]
        dot << "digraph structs {\n  newrank=true;\n  dpi=72;\n";
        dot << "  subgraph cluster1 { style=invis;\n";
        dot << "  node [shape=record, height=0.8, width=0.5, fontsize=50, penwidth=4, fillcolor=white, style=filled];\n";
        dot << "  edge [arrowhead=normal, minlen=1, penwidth=4, color=blue];\n";
        const int max_levels = abs(e.getLevel());
        nodes_per_lvl.resize(num_series);
        for (size_t i=0; i<num_series; i++)
            nodes_per_lvl[i].resize(max_levels + 1);

        visit(e.getNode());

        dot << "}\n";

        if (write_level_labels) {
            // write places/levels in a separate cluster
            dot << "  subgraph cluster2 { style=invis;\n  node [shape=none, fontsize=60, margin=\"0.5,0.1\"];\n";
            for (int lvl=1; lvl<=max_levels; lvl++) {
                const char *level_name = forest->getDomain()->getVar(lvl)->getName();
                for (size_t i=0; i<num_series; i++) {
                    dot << "  LVL"<<lvl<<"_"<<i<<" [label=\""
                        << level_name<<(i==1?"\\\'":"")<<"\"];\n";
                    nodes_per_lvl[i][lvl] += " LVL";
                    nodes_per_lvl[i][lvl] += std::to_string(lvl);
                    nodes_per_lvl[i][lvl] += "_";
                    nodes_per_lvl[i][lvl] += std::to_string(i);
                }
            }
            if (write_terminals)
                dot << "  LVL_terms [label=\"\"];\n";
            dot << "}\n";
        }

        // write node rankings
        for (int lvl = max_levels; lvl > 0; lvl--) {
            for (size_t i=0; i<num_series; i++) {
                dot << (lvl==max_levels && i==0 ? "" : " -> ") 
                    << "{rank=same " << nodes_per_lvl[i][lvl] << "}";
            }
        }
        if (write_terminals)
            dot << " -> {rank=same  nT LVL_terms}";
        dot << " [style=invis]\n";
        dot << "}\n";
    }
};

/////////////////////////////////////////////////////////////////////////////////////////

void write_dd_as_dot(const MEDDLY::dd_edge& e, 
                     const char* dot_name, bool level_labels, 
                     bool write_terminals)
{
    dot_of_DD d;
    d.write_level_labels = level_labels;
    d.write_terminals    = write_terminals;
    d.dot.open(dot_name);
    d.start_visit(e);
}

/////////////////////////////////////////////////////////////////////////////////////////

int dot_to_pdf(const char *dot_fname, const char *pdf_fname)
{
    char buffer[256];
    snprintf(buffer, sizeof(buffer), "dot \"%s\" -Tpdf -o \"%s\"", dot_fname, pdf_fname);
    return system(buffer);
    // const char* const args[] = { "dot", dot_fname, "-Tpdf", "-o", pdf_fname, NULL };
    // system();
}

/////////////////////////////////////////////////////////////////////////////////////////

void write_dd_as_pdf(const MEDDLY::dd_edge& e, 
                     const char* base_name, bool level_labels, 
                     bool write_terminals)
{
    std::string dot_name(base_name), pdf_name(base_name);
    dot_name += ".dot";
    pdf_name += ".pdf";

    write_dd_as_dot(e, dot_name.c_str(), level_labels, write_terminals);

    dot_to_pdf(dot_name.c_str(), pdf_name.c_str());
}

/////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>

#include "sym_pottier.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////

// Save an MDD as an explicit matrix
void save_mat_explict(const std::string& base_fname, const char* ext,
                      const variable_order& vorder, MEDDLY::dd_edge dd) 
{
    std::string mat_fname = base_fname + ext;
    cout << "Saving " << mat_fname << " ..." << endl;
    ofstream ofs(mat_fname);
    ofs << print_mdd(dd, vorder, false, true) <<endl;
    ofs.close();
}

/////////////////////////////////////////////////////////////////////////////////////////

// Serialize an MDD using the MEDDLY format
void save_mat_dd(const std::string& base_fname, const char* ext,
                 MEDDLY::forest* forestMDD, MEDDLY::dd_edge dd) 
{
    std::string mat_fname = base_fname + ext;
    cout << "Saving " << mat_fname << " ..." << endl;
    ofstream ofs(mat_fname);
    MEDDLY::ostream_output mofs(ofs);

    MEDDLY::dd_edge list[1];
    list[0] = dd;
    forestMDD->writeEdges(mofs, list, 1);
    mofs.flush();
}

/////////////////////////////////////////////////////////////////////////////////////////

#if (defined __unix__) || (defined __APPLE__)

#include <sys/resource.h>
void print_cpu_rusage() {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    cout << "CPU usage: " << r.ru_utime.tv_sec << "." 
         << setw(6) << setfill('0') << r.ru_utime.tv_usec << setfill(' ') << endl;
}

#else

void print_cpu_rusage() {
    cout << "CPU rusage unsupported.";
}

#endif

/////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////

// void test() {
//     const std::vector<std::vector<int>> mA = {
//         { 0, 3, 1, 0, 0},
//         { 2, 0, 0, 1, 2},
//         { 1,-1,-3,-1,-1},
//         { 0, 2, 0, 2, 0},
//         { 0, 2, 4, 2, 0},
//         { 0, 2, 3, 2, 0},
//     };
//     const std::vector<std::vector<int>> mB = {
//         { 3,-2,-2, 1, 4},
//         {-5, 1, 1, 2,-3},
//         { 0, 1, 3, 0, 0},
//     };
//     const size_t num_levels = mA[0].size();
//     variable_order vorder(num_levels), pivot_order(num_levels);

//     meddly_context ctx(num_levels, vorder, pivot_order);
//     size_t meddly_cache_size = 1000000;
//     ctx.initialize(meddly_cache_size);
//     MEDDLY::dd_edge A = mdd_from_vectors(mA, ctx.forestMDD, false);
//     MEDDLY::dd_edge B = mdd_from_vectors(mB, ctx.forestMDD, false);

//     cout << "A:\n" << print_mdd(A, vorder) << endl;
//     // cout << "B:\n" << print_mdd(B, vorder) << endl;

//     MEDDLY::dd_edge U(ctx.forestMDD), D(ctx.forestMDD);
//     VCANON2->computeDDEdge(A, 2, U, D);
//     cout << "U:\n" << print_mdd(U, vorder) << endl;
//     cout << "D:\n" << print_mdd(D, vorder) << endl;

//     // MEDDLY::dd_edge AB(ctx.forestMDD);
//     // s_vectors* svop = S_VECTORS_OPS->get_op(1, false, false, SVS_UNDECIDED); // lvl, isConf, invertB, svs
//     // svop->computeDDEdge(A, B, AB, false);
//     // cout << "A + B:\n" << print_mdd(AB, vorder) << endl;

//     // MEDDLY::dd_edge A2(ctx.forestMDD);
//     // SIGN_CANON_OPS->get_op(true)->computeDDEdge(A, A2, false);
//     // cout << "A:\n" << print_mdd(A2, vorder) << endl;

//     // MEDDLY::dd_edge AB = sym_s_vectors_extreme_rays(A, B, 3, vorder);
//     // cout << "A+B:\n" << print_mdd(AB, vorder) << endl;

//     // divisors_finder_mdd_op::lut_key lk;
//     // DIV_FINDER_MDD->compute(A, lk);
//     // const std::vector<int>& divisors = DIV_FINDER_MDD->look_up(lk);

//     // cout << print_mdd(A, vorder) << endl << "Divisors: ";

//     // // divisors appear in descending order
//     // for (int div : divisors)
//     //     cout << div << endl;
//     // cout << endl << endl;
//     exit(0);
// }

/////////////////////////////////////////////////////////////////////////////////////////

// void test2() {
//     std::vector<std::vector<int>> A, kerA;
//     read_mat("test/k44bad.mat", A);
//     cout << "A:" << endl; print_mat(A); cout << endl;

//     // integral_kernel(A, kerA, false);
//     // cout << "kerA 1:" << endl; print_mat(kerA); cout << endl;

//     std::vector<std::vector<int>> H, U, UA;
//     std::vector<size_t> leading_cols;
//     hermite_normal_form(transpose(A), H, U, &leading_cols, true);
//     cout << leading_cols.size() << endl;

//     // kerA = integral_kernel(A);
//     // cout << "kerA 2:" << endl; print_mat(kerA); cout << endl;


//     // std::vector<std::vector<int>> ief = integral_echelon_form(A, nullptr);
//     // cout << "ief(A):" << endl; print_mat(ief); cout << endl;

//     exit(0);
// }


// void test3() {
//     MEDDLY::initializer_list* L = MEDDLY::defaultInitializerList(0);
//     MEDDLY::initialize(L);

//     // Initialize domain
//     MEDDLY::domain *dom;
//     size_t num_levels = 5;
//     variable_order vorder(num_levels);
//     vector<int> domainBnd(num_levels);
//     for (int i=0; i<num_levels; i++)
//         domainBnd[i] = 1;
//     dom = MEDDLY::createDomainBottomUp(domainBnd.data(), num_levels);
//     for (int i=0; i<num_levels; i++) {
//         char buffer[64];
//         snprintf(buffer, sizeof(buffer), "x%lu", vorder.lvl2var(i)+1);
//         dom->useVar(i + 1)->setName(strdup(buffer));
//     }

//     // Initialize the MDD forest
//     MEDDLY::policies mdd_fp(false);
//     mdd_fp.setQuasiReduced();
//     mdd_fp.setSparseStorage();
//     MEDDLY::forest *forestMDD;
//     forestMDD = dom->createForest(false, MEDDLY::range_type::BOOLEAN, 
//                                   MEDDLY::edge_labeling::MULTI_TERMINAL, mdd_fp);

//     const std::vector<std::vector<int>> mA = {
//         { 0, 0, 0,-1, 1},
//     };
//     const std::vector<std::vector<int>> mB = {
//         { 0, 0, -1, 0, 1},
//     };

//     // meddly_context ctx(num_levels, vorder, pivot_order);
//     MEDDLY::dd_edge A = mdd_from_vectors(mA, forestMDD, false);
//     MEDDLY::dd_edge B = mdd_from_vectors(mB, forestMDD, false);

//     cout << "A:\n" << print_mdd(A, vorder) << endl;
//     cout << "B:\n" << print_mdd(B, vorder) << endl;    

//     MEDDLY::dd_edge AuB(forestMDD);
//     MEDDLY::apply(MEDDLY::UNION, A, B, AuB);

//     cout << "A U B:\n" << print_mdd(AuB, vorder) << endl;    

//     exit(0);
// }

/////////////////////////////////////////////////////////////////////////////////////////
// Unit tests:

void unit_test_vcanon();
void unit_test_reduce();
void unit_test_sign_canon();
void unit_test_minimal_supports();

void do_unit_tests(int argc, char** argv) {
    int ii=1;
    while (ii < argc) {
        if (0==strcmp(argv[ii], "-unit-test-vcanon")) { 
            unit_test_vcanon(); 
        }
        else if (0==strcmp(argv[ii], "-unit-test-reduce")) { 
            unit_test_reduce(); 
        }
        else if (0==strcmp(argv[ii], "-unit-test-sign-canon")) { 
            unit_test_sign_canon(); 
        }
        else if (0==strcmp(argv[ii], "-unit-test-minimal-supports")) { 
            unit_test_minimal_supports(); 
        }
        ++ii;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    do_unit_tests(argc, argv);
    // test();
    // test2();
    // test3();
    const char* basename = nullptr;
    bool transpose_input = false;
    bool do_check_output = false;
    selected_varorder svo = selected_varorder::SLOAN;
    bool compute_Zgenerators = true;
    bool show_help = false;
    bool save_dd = false;
    bool output_mat_explicit = false;
    bool output_mat_dd = false;
    bool do_pivoting = true;
    bool print_rusage = false;
    pottier_params_t pparams;
    perf_counters_t perf_counters;
    bool use_old_basisZgen = false;
    size_t num_iters_pivot_heur = 2;
    size_t meddly_cache_size = 1000000;

    // Read command-line arguments
    int ii=1;
    while (ii < argc) {
        if (0==strcmp(argv[ii], "-v")) {
            pparams.verbose = true;
        }
        else if (0==strcmp(argv[ii], "-vb")) {
            pparams.verbose = pparams.verbose_for_basis = true;
        }
        else if (0==strcmp(argv[ii], "-vv")) {
            pparams.verbose = pparams.very_verbose = pparams.verbose_for_basis = true;
        }
        else if (0==strcmp(argv[ii], "-t")) {
            transpose_input = true;
        }
        else if (0==strcmp(argv[ii], "-k")) {
            compute_Zgenerators = false;
        }
        else if (0==strcmp(argv[ii], "-yl")) {
            pparams.by_levels = true;
            pparams.normalize_by_levels = false;
        }
        else if (0==strcmp(argv[ii], "-nl")) {
            pparams.by_levels = false;
        }
        // else if (0==strcmp(argv[ii], "-z")) { // TODO: disable
        //     pparams.by_levels = true;
        //     pparams.normalize_by_levels = true;
        // }
        else if (0==strcmp(argv[ii], "-s")) { 
            pparams.by_levels = true;
            pparams.dynamic_svectors = true;
        }
        else if (0==strcmp(argv[ii], "-ye")) {
            pparams.by_generators = true;
        }
        else if (0==strcmp(argv[ii], "-ne")) {
            pparams.by_generators = false;
        }
        else if (0==strcmp(argv[ii], "-p")) {
            pparams.perf = &perf_counters;
        }
        else if (0==strcmp(argv[ii], "-g")) {
            pparams.target = compute_target::GRAVER_BASIS;
        }
        else if (0==strcmp(argv[ii], "-x")) {
            pparams.target = compute_target::EXTREME_RAYS;
        }
        else if (0==strcmp(argv[ii], "-cs") && ii+1 < argc) {
            meddly_cache_size = stoul(argv[++ii]);
        }
        else if (0==strcmp(argv[ii], "-yd")) {
            pparams.by_levels = true;
            pparams.by_degree = true;
        }
        else if (0==strcmp(argv[ii], "-nd")) {
            pparams.by_degree = false;
        }
        else if (0==strcmp(argv[ii], "-c")) {
            do_check_output = true;
        }
        else if (0==strcmp(argv[ii], "-no")) {
            svo = selected_varorder::NONE;
        }
        else if (0==strcmp(argv[ii], "-hn")) { // TODO: remove
            use_old_basisZgen = false;
        }
        else if (0==strcmp(argv[ii], "-ige")) { // TODO: remove
            use_old_basisZgen = true;
        }
        else if (0==strcmp(argv[ii], "-oldp")) { // TODO: remove option
            num_iters_pivot_heur = 1;
        }
#ifdef HAS_BOOST_CPP
        else if (0==strcmp(argv[ii], "-sloan2")) {
            svo = selected_varorder::BOOST_SLOAN;
        }
#endif // HAS_BOOST_CPP
        else if (0==strcmp(argv[ii], "-sloan")) {
            svo = selected_varorder::SLOAN;
        }
        else if (0==strcmp(argv[ii], "-fast")) {
            svo = selected_varorder::FAST;
        }
        else if (0==strcmp(argv[ii], "-po")) {
            svo = selected_varorder::PIVOTING;
        }
        else if (0==strcmp(argv[ii], "-dd")) {
            save_dd = true;
        }
        else if (0==strcmp(argv[ii], "-o")) {
            output_mat_explicit = true;
        }
        else if (0==strcmp(argv[ii], "-od")) {
            output_mat_dd = true;
        }
        else if (0==strcmp(argv[ii], "-np")) {
            do_pivoting = false;
        }
        else if (0==strcmp(argv[ii], "-u")) {
            print_rusage = true;
        }
        else if (0==strcmp(argv[ii], "-h")) {
            show_help = true;
        }
        else if (basename==nullptr) {
            basename = argv[ii];
        }
        else {
            cerr << "Unknown option: " << argv[ii] << endl;
            return 1;
        }
        ii++;
    }
    if (basename == nullptr)
        show_help = true;

    if (show_help) {
        cout << ("sym_hilbert: compute the Hilbert basis using Decision Diagrams.\n"
                 "  Usage: sym_hilbert <options> <basename>\n"
                 "  Input matrix should be called <basename>.mat\n\n"
                 "Options: -v       Verbose.\n"
                 "         -vv      Very Verbose.\n"
                 "         -p       Performance counters.\n"
                 "         -u       Print CPU rusage.\n"
                 "         -t       Transpose input matrix.\n"
                 "         -k       Input matrix is already the lattice generating set.\n"
                 "         -g       Compute the Graver basis instead of the Hilbert basis.\n"
                 "         -x       Compute the set of extremal rays of the cone.\n"
                 "         -cs <x>  Set Meddly cache size as <x>.\n"
                 "Method:  -yl      Computate by levels. [default]\n"
                 "         -nl      Disable computation by levels.\n"
                 "         -z       Compute by levels and normalize by levels.\n"
                 "         -yd      Evaluate candidates in graded order (by levels) [default].\n"
                 "         -nd      Do not evaluate candidates in graded order (by levels).\n"
                 "         -s       Generate candidates in graded order (by levels).\n"
                 "         -ye      Compute by generators.\n"
                 "         -ne      Disable computation by generators. [default]\n"
                 "         -np      Disable level pivoting.\n"
                 "Output:  -dd      Save decision diagrams as PDFs.\n"
                 "         -o       Write results as explicit matrices.\n"
                 "         -od      Write results as serialized MDDs (MEDDLY format).\n"
                 "Reorder: -no      Do not reorder the problem variables.\n"
#ifdef HAS_BOOST_CPP
                 "         -sloan2  Use Sloan algorithm (boost c++ version) to reorder variables.\n"
#endif // HAS_BOOST_CPP
                 "         -sloan   Use Sloan algorithm to reorder variables (default).\n"
                 "         -fast    Use fast algorithm to reorder variables.\n"
                 "         -po      Use pivoting order for both MDD levels and pivoting.\n"
                 ) << endl;
    }
    if (basename == nullptr) {
        cerr << "Missing input filename." << endl;
        return 1;
    }
    if (pparams.target == compute_target::EXTREME_RAYS) {
        pparams.by_levels = true; // extreme rays can only be computed by variables
    }
    if (pparams.target == compute_target::GRAVER_BASIS) {
        pparams.dynamic_svectors = false; // degree is defined only for Hilbert/rays
    }
    if (pparams.target != compute_target::HILBERT_BASIS || pparams.by_generators) {
        pparams.by_degree = false; // only for Hilbert bases
    }
    if (pparams.by_generators) {
        pparams.normalize_by_levels = false;
    }
    if (!pparams.by_levels) {
        pparams.by_degree = false;
    }

    const std::string base_fname(basename);

    //----------------------------------------------------
    // Setup of the linear eq. system
    //----------------------------------------------------

    // Read input matrix
    std::string mat_fname = base_fname + (compute_Zgenerators ? ".mat" : ".lat");
    std::vector<std::vector<int>> A;
    if (!read_mat(mat_fname.c_str(), A)) {
        cerr << "Cannot read input matrix file: " << mat_fname << endl;
        return 1;
    }
    if (A.size() == 0) {
        cerr << "Input matrix is empty." << endl;
        return 0;
    }
    if (pparams.verbose_show_mat(A)) {
        cout << "Input matrix:" << endl; print_mat(A); cout << endl;
    }

    // Transpose
    if (transpose_input) {
        A = transpose(A);
    }

    const size_t num_variables = A.front().size();

    // {
    //     std::vector<std::vector<int>> ker1A, ker2A;

    //     integral_kernel_old(A, ker1A, true); 

    //     cout << "\n\n\n\nNEW METHOD:\n";
    //     ker2A = integral_kernel_Zgens(A);

    //     cout << endl << endl << endl << endl << endl;
    //     cout << "Z-basis for the integral kernel of A (OLD):" << endl; 
    //     print_mat(ker1A); cout << endl;
    //     cout << "Z-basis for the integral kernel of A (NEW):" << endl; 
    //     print_mat(ker2A); cout << endl;

    //     if (ker1A == ker2A) {
    //         cout << "OK" << endl;
    //         exit(0);
    //     }
    //     else {
    //         cout << "NOT SAME" << endl;
    //         exit(2);
    //     }
    // }

    // Compute the Z-generators of the lattice
    std::vector<std::vector<int>> lattice_Zgenerators;
    if (compute_Zgenerators) {
        if (use_old_basisZgen) {
            // TODO: remove
            integral_kernel_old(A, lattice_Zgenerators, pparams.verbose_for_basis); 
        }
        else {
            lattice_Zgenerators = integral_kernel_Zgens(A, pparams.verbose_for_basis);
        }

        if (pparams.verbose_show_mat(lattice_Zgenerators)) {
            cout << "Z-basis for the integral kernel of A:" << endl; print_mat(lattice_Zgenerators); cout << endl;
        }
    }
    else {
        lattice_Zgenerators = A;
    }
    if (lattice_Zgenerators.size() == 0) {
        cerr << "Z-basis is empty." << endl;
        // return 0;
    }

    //----------------------------------------------------
    // Variable ordering to minimize DD size
    //----------------------------------------------------

    // Select the matrix that will be used to compute the variable order.
    // In principle it should be the constraints matrix A, 
    // but if it is not available, we will use lattice_Zgenerators.
    const std::vector<std::vector<int>> *p_reorder_mat; // problem constraints
    p_reorder_mat = compute_Zgenerators ? &A : &lattice_Zgenerators;

    // Find a reasonable variable order for this problem
    variable_order vorder(num_variables);

    switch (svo) {
        case selected_varorder::NONE:
            break;
#ifdef HAS_BOOST_CPP
        case selected_varorder::BOOST_SLOAN:
            boost_sloan_varorder(*p_reorder_mat, vorder);
            break;
#endif // HAS_BOOST_CPP
        case selected_varorder::SLOAN:
            sloan_varorder(*p_reorder_mat, vorder);
            break;
        case selected_varorder::FAST:
            fast_varorder(*p_reorder_mat, vorder);
            break;
        case selected_varorder::PIVOTING:
            pivot_order_from_matrix_iter(vorder, lattice_Zgenerators, pparams.target==compute_target::GRAVER_BASIS, 
                                         num_iters_pivot_heur);
            vorder.invert();
            do_pivoting = false;
            break;
    }

    if (svo != selected_varorder::NONE) {
        if (pparams.verbose) {
            cout << "MDD variable order:\n";
            vorder.print();
        }
        if (pparams.very_verbose)
            cout << "Initial iRank: " << irank(*p_reorder_mat) << endl;
        if (pparams.very_verbose) {
            std::vector<std::vector<int>> reorderedA = reorder_matrix(*p_reorder_mat, vorder);
            row_footprint_form(reorderedA);
            cout << "iRank after reordering: " << irank(reorderedA) << endl;
            cout << "Reordered Input matrix:" << endl; print_mat(reorderedA, true); cout << endl;
        }

        lattice_Zgenerators = reorder_matrix(lattice_Zgenerators, vorder);
    }

    // Pivot ordering
    variable_order pivot_order(num_variables);
    if (do_pivoting) {
        // pivot_order_from_matrix(pivot_order, lattice_Zgenerators, pparams.target==compute_target::GRAVER_BASIS);
        pivot_order_from_matrix_iter(pivot_order, lattice_Zgenerators, pparams.target==compute_target::GRAVER_BASIS, 
                                     num_iters_pivot_heur);
    }

    if (pparams.verbose) {
        cout << "Pivot order:\n";
        pivot_order.print();

        // for (size_t var=0; var<m; var++) {
        //     size_t lambda = pivot_order.var2lvl(var) + 1;
        //     for (size_t l=1; l<=m; l++) {
        //         if (pivot_order.is_same_as_lambda(lambda, l)) cout << "=";
        //         if (pivot_order.is_above_lambda(lambda, l)) cout << ".";
        //         if (pivot_order.is_below_lambda(lambda, l)) cout << "*";
        //     }
        //     cout << endl;
        // }
        // cout << endl;
    }

    // canonicalize_by_order(lattice_Zgenerators, pivot_order);

        // auto rrffKerA = lattice_Zgenerators;
        // row_footprint_form(rrffKerA);
        // cout << "RRFF(Kernel(A)):" << endl; print_mat(rrffKerA, true); cout << endl;
    // exit(0);

    // Initialize MEDDLY and create an MDD forest for the computations
    meddly_context ctx(num_variables, vorder, pivot_order);
    ctx.initialize(meddly_cache_size);

    //----------------------------------------------------
    // Symbolic computation of the Graver basis
    //----------------------------------------------------
    sep(pparams);

    // Get the Graver basis using the Pottier algorithm
    MEDDLY::dd_edge G = sym_pottier_bygen(ctx, pparams, lattice_Zgenerators);

    sep(pparams);

    //----------------------------------------------------
    // Show and save results
    //----------------------------------------------------
    struct output_params {
        const char *what, *mat_ext, *dd_ext, *pdf_suffix;
    } op;
    switch (pparams.target) {
        case compute_target::GRAVER_BASIS:
            op = output_params{ 
                .what="Graver basis", .mat_ext = ".grav", 
                .dd_ext=".gravdd", .pdf_suffix="-graver" };
            break;
        case compute_target::HILBERT_BASIS:
            op = output_params{ 
                .what="Hilbert basis", .mat_ext = ".hilb", 
                .dd_ext=".hilbdd", .pdf_suffix="-hilbert" };
            break;
        case compute_target::EXTREME_RAYS:
            op = output_params{ 
                .what="Extreme rays set", .mat_ext = ".xray", 
                .dd_ext=".xraydd", .pdf_suffix="-exrays" };
            break;
    }
    double card_G = G.getCardinality();
    cout << "\n"<<op.what<<" has "<< card_G << " entries.\n" << endl;
    if (pparams.very_verbose || (pparams.verbose && card_G<50)) {
        cout << print_mdd(G, vorder) << endl;
    }
    if (save_dd) {
        std::string dd_fname = base_fname + op.pdf_suffix;
        cout << "Saving " << dd_fname << ".pdf" << endl;
        write_dd_as_pdf(G, dd_fname.c_str(), true, true);
    }
    if (output_mat_explicit) {
        save_mat_explict(base_fname, op.mat_ext, vorder, G);
    }
    if (output_mat_dd) {
        save_mat_dd(base_fname, op.dd_ext, ctx.forestMDD, G);
    }
    if (pparams.perf) {
        cout << "Total S-Vectors processed: "<<pparams.perf->counter_C<<endl;
    }

    if (print_rusage)
        print_cpu_rusage();

    //----------------------------------------------------
    // Check correctness of generated solutions.
    //----------------------------------------------------
    if (do_check_output) {
        std::vector<std::vector<int>> Hcheck; // results to check
        bool can_check = false, canonicalize_half_basis = false;
        switch (pparams.target) {
            case compute_target::GRAVER_BASIS: {
                std::string gra_fname = base_fname + ".gra";
                if (!read_mat(gra_fname.c_str(), Hcheck)) {
                    cerr << "Cannot read Graver matrix file: " << gra_fname << endl;
                }
                else {
                    can_check = true;
                    canonicalize_half_basis = true;
                }
            }
            break;
            case compute_target::HILBERT_BASIS: {
                std::string hil_fname = base_fname + ".hil";
                if (!read_mat(hil_fname.c_str(), Hcheck)) {
                    cerr << "Cannot read Hilbert matrix file: " << hil_fname << endl;
                }
                else can_check = true;
            }
            break;
            case compute_target::EXTREME_RAYS: {
                std::string ray_fname = base_fname + ".ray";
                if (!read_mat(ray_fname.c_str(), Hcheck)) {
                // std::string ray_fname = base_fname + ".pin";
                // if (!read_pin(ray_fname.c_str(), m, Hcheck)) {
                    cerr << "Cannot read Extreme Rays matrix file: " << ray_fname << endl;
                }
                else can_check = true;
            }
            break;
        }

        if (can_check) {
            Hcheck = reorder_matrix(Hcheck, vorder);

            MEDDLY::dd_edge Hcheckdd = mdd_from_vectors(Hcheck, ctx.forestMDD, false);
            if (canonicalize_half_basis) {
                SIGN_CANON->computeDDEdge(Hcheckdd, Hcheckdd, false);
            }

            if (Hcheckdd == G) {
                cout << "[OK]" << endl;
            }
            else {
                cout << "[FAILED]" << endl;

                if (pparams.verbose) {
                cout << "Expected size: " << Hcheckdd.getCardinality() << endl;
                    MEDDLY::dd_edge diff(ctx.forestMDD);
                    // MEDDLY::apply(MEDDLY::DIFFERENCE, G, Hcheckdd, diff);
                    diff = sym_difference(G, Hcheckdd);
                    cout << "G \\ Expected:\n" << print_mdd(diff, ctx.vorder) << endl;
                    // MEDDLY::apply(MEDDLY::DIFFERENCE, Hcheckdd, G, diff);
                    diff = sym_difference(Hcheckdd, G);
                    cout << "Expected \\ G:\n" << print_mdd(diff, ctx.vorder) << endl;
                }

                return 2;
            }
        }
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////





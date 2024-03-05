/////////////////////////////////////////////////////////////////////////////////////////
#ifndef __VARIABLE_ORDER_H__
#define __VARIABLE_ORDER_H__
/////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
// Variable reordering algorithms
/////////////////////////////////////////////////////////////////////////////////////////

class variable_order {
    std::vector<size_t> var_to_level, level_to_var;
public:
    variable_order(size_t m);

    inline size_t var2lvl(size_t var) const { return var_to_level[var]; }
    inline size_t lvl2var(size_t lvl) const { return level_to_var[lvl]; }

    inline void bind_var2lvl(size_t var, size_t lvl) {
        var_to_level[var] = lvl; level_to_var[lvl] = var;
    }

    void check_order() const;

    void print() const;

    inline void invert() { var_to_level.swap(level_to_var); }

    // Pivoting for meddly
    // same as lambda means that lvl == lambda in pivot order
    inline bool is_same_as_lambda(size_t lambda, size_t lvl) const 
    {   assert(lvl-1 < level_to_var.size() && lambda-1 < level_to_var.size()); 
        return lvl2var(lvl-1)==lvl2var(lambda-1);
    }
    // above lambda means that lvl has not been encountered yet in the pivot order
    inline bool is_above_lambda(size_t lambda, size_t lvl) const 
    {   assert(lvl-1 < level_to_var.size() && lambda-1 < level_to_var.size()); 
        return lvl2var(lvl-1) > lvl2var(lambda-1);
    }
    // below lambda means that lvl has already been processed in the pivot order
    inline bool is_below_lambda(size_t lambda, size_t lvl) const 
    {   assert(lvl-1 < level_to_var.size() && lambda-1 < level_to_var.size()); 
        return lvl2var(lvl-1) < lvl2var(lambda-1);
    }
};

// reorder columns of a matrix
std::vector<std::vector<int>>
reorder_matrix(const std::vector<std::vector<int>>& A,
               const variable_order &vorder);

// // ensure that first non-zero entry in that specific order is positive
// void canonicalize_by_order(std::vector<std::vector<int>>& mat,
//                            const variable_order &vorder);

#ifdef HAS_BOOST_CPP
void boost_sloan_varorder(const std::vector<std::vector<int>>& A,
                          variable_order &out_order);
#endif

void fast_varorder(const std::vector<std::vector<int>>& A,
                   variable_order &out_order);

void sloan_varorder(const std::vector<std::vector<int>>& A,
                    variable_order &out_order);

enum class selected_varorder {
    NONE, BOOST_SLOAN, SLOAN, FAST, PIVOTING
};

// pivot ordering for Pottier-by-level
// void pivot_order_from_matrix(variable_order& pivots,
//                              const std::vector<std::vector<int>>& A,
//                              bool optimize_graver);

void pivot_order_from_matrix_iter(variable_order& pivots,
                                  const std::vector<std::vector<int>>& A,
                                  const bool optimize_graver,
                                  const size_t num_iters,
                                  const std::vector<size_t>& fixed_vars);

/////////////////////////////////////////////////////////////////////////////////////////
#endif // __VARIABLE_ORDER_H__
/////////////////////////////////////////////////////////////////////////////////////////

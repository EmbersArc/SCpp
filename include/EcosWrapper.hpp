#pragma once

#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <sstream> 
#include <experimental/optional> 
#include <utility>
#include <functional>
#include "OptimizationProblem.hpp"

#include "ecos.h"
using std::vector;
using std::string;
using std::map;
using std::experimental::optional;
using std::pair;


void sparse_DOK_to_CCS(
    const map< pair<idxint, idxint>, optimization_problem::Parameter > &sparse_DOK, 
    vector<optimization_problem::Parameter> &data_CCS, 
    vector<idxint> &columns_CCS, 
    vector<idxint> &rows_CCS,
    size_t n_columns);



class EcosWrapper {



    /* ECOS problem parameters */
    idxint                                  ecos_n_variables;
    idxint                                  ecos_n_constraint_rows;
    idxint                                  ecos_n_equalities;
    idxint                                  ecos_n_positive_constraints;
    idxint                                  ecos_n_cone_constraints;
    vector<idxint>                          ecos_cone_constraint_dimensions;
    idxint                                  ecos_n_exponential_cones;
    vector<optimization_problem::Parameter> ecos_G_data_CCS;
    vector<idxint>                          ecos_G_columns_CCS;
    vector<idxint>                          ecos_G_rows_CCS;
    vector<optimization_problem::Parameter> ecos_A_data_CCS;
    vector<idxint>                          ecos_A_columns_CCS;
    vector<idxint>                          ecos_A_rows_CCS;
    vector<optimization_problem::Parameter> ecos_cost_function_weights;
    vector<optimization_problem::Parameter> ecos_h;
    vector<optimization_problem::Parameter> ecos_b;

    /* ECOS result */
    vector<double> ecos_solution_vector;
    idxint ecos_exitflag;


public:
    optimization_problem::SecondOrderConeProgram socp;

    explicit EcosWrapper(optimization_problem::SecondOrderConeProgram _socp);
    
    //void compile_problem_structure();
    void solve_problem();

    double get_solution_value(size_t problem_index) {
        return ecos_solution_vector[problem_index];
    }

    double get_solution_value(const string &name, const vector<size_t> &indices) {
        return ecos_solution_vector[socp.get_variable(name, indices).problem_index];
    }
};
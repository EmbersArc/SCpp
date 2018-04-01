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
    size_t n_variables = 0;

    /* Set of named tensor variables in the optimization problem */
    map<string, vector<size_t>> tensor_variable_dimensions;
    map<string, vector<size_t>> tensor_variable_indices;

    /* High-level optimization problem description */
    vector<optimization_problem::SecondOrderConeConstraint> secondOrderConeConstraints;
    vector<optimization_problem::PostiveConstraint> postiveConstraints;
    vector<optimization_problem::EqualityConstraint> equalityConstraints;
    optimization_problem::AffineExpression costFunction = optimization_problem::Parameter(0.0);

    /* ECOS problem parameters */
    idxint         ecos_n_variables;
    idxint         ecos_n_constraint_rows;
    idxint         ecos_n_equalities;
    idxint         ecos_n_positive_constraints;
    idxint         ecos_n_cone_constraints;
    vector<idxint> ecos_cone_constraint_dimensions;
    idxint         ecos_n_exponential_cones;
    vector<optimization_problem::Parameter> ecos_G_data_CCS;
    vector<idxint> ecos_G_columns_CCS;
    vector<idxint> ecos_G_rows_CCS;
    vector<optimization_problem::Parameter> ecos_A_data_CCS;
    vector<idxint> ecos_A_columns_CCS;
    vector<idxint> ecos_A_rows_CCS;
    vector<optimization_problem::Parameter> ecos_cost_function_weights;
    vector<optimization_problem::Parameter> ecos_h;
    vector<optimization_problem::Parameter> ecos_b;

    /* ECOS result */
    vector<double> ecos_solution_vector;
    idxint ecos_exitflag;

    size_t allocate_variable_index();

public:
    void create_tensor_variable(const string &name, const vector<size_t> &dimensions);
    size_t get_tensor_variable_index(const string &name, const vector<size_t> &indices);
    optimization_problem::Variable get_variable(const string &name, const vector<size_t> &indices);

    void add_constraint(optimization_problem::SecondOrderConeConstraint c);
    void add_constraint(optimization_problem::PostiveConstraint c);
    void add_constraint(optimization_problem::EqualityConstraint c);
    void add_minimization_term(optimization_problem::AffineExpression c);
    void compile_problem_structure();
    void solve_problem();

    double get_solution_value(size_t problem_index) {
        return ecos_solution_vector[problem_index];
    }

    double get_solution_value(const string &name, const vector<size_t> &indices) {
        return ecos_solution_vector[get_variable(name, indices).problem_index];
    }

    void print_problem(std::ostream &out);
};
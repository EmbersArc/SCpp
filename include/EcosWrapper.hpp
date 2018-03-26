#pragma once

#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <sstream> 
#include <experimental/optional> 
#include "ecos.h"
using std::vector;
using std::string;
using std::map;
using std::experimental::optional;

namespace optimization_problem {
    struct Variable { // represents an optimization variable x_i
        string name;
        vector<size_t> tensor_indices;
        size_t problem_index;
        string print() const;
    };

    struct AffineTerm;

    struct AffineExpression { // represents a term like (p_1*x_1 + p_2*x_2 + b)
        vector<AffineTerm> terms;
        string print() const;
    };

    class Parameter { // represents a parameter p_i in the opt-problem that can be changed between problem evaluations.
        const double* dynamic_value_ptr = NULL;
        double const_value = 0;
        bool is_const;
    public:
        Parameter(const double* dynamic_value_ptr):dynamic_value_ptr(dynamic_value_ptr),is_const(false) {
            if(dynamic_value_ptr == NULL) throw std::runtime_error("Null Pointer Error: Parameter(NULL)");
        }
        Parameter(double const_value):const_value(const_value),is_const(true){}
        Parameter():const_value(0),is_const(true){}
        double get_value() const {
            if(is_const) return const_value;
            return *dynamic_value_ptr;
        }
        string print() const;
        operator AffineTerm();
        operator AffineExpression();
    };

    struct AffineTerm { // represents a linear term (p_1*x_1) or constant term (p_1)
        Parameter parameter = Parameter(0.0);
        optional<Variable> variable;
        string print() const;
        operator AffineExpression();
    };

    struct Norm2 { // represents a term like norm2([p_1*x_1 + p_2*x_2 + b_1,   p_3*x_3 + p_4*x_4 + b_2 ])
        vector<AffineExpression> arguments;
        string print() const;
    };

    // represents a constraint like 
    //      norm2([p_1*x_1 + p_2*x_2 + b_1,   p_3*x_3 + p_4*x_4 + b_2 ])
    //        <= p_5*x_5 + p_6*x_6 + b_3
    struct SecondOrderConeConstraint {
        Norm2 lhs;
        AffineExpression rhs;
        string print() const;
    };


    // represents a constraint like 
    //     p_1*x_1 + p_2*x_2 + b >= 0
    struct PostiveConstraint
    {
        AffineExpression lhs;
        string print() const;
    };

    // represents a constraint like 
    //     p_1*x_1 + p_2*x_2 + b == 0
    struct EqualityConstraint
    {
        AffineExpression lhs;
        string print() const;
    };

    AffineTerm operator*(const Parameter &parameter, const Variable &variable);
    AffineTerm operator*(const double &const_parameter, const Variable &variable);
    AffineExpression operator+(const AffineExpression &lhs, const AffineExpression &rhs);
    AffineExpression operator+(const AffineExpression &lhs, const double &rhs);
    Norm2 norm2(const vector<AffineExpression> &affineExpressions);
    SecondOrderConeConstraint operator<=(const Norm2 &lhs, const AffineExpression &rhs);
    SecondOrderConeConstraint operator<=(const Norm2 &lhs, const double &rhs);
    PostiveConstraint operator>=(const AffineExpression &lhs, const double &zero);
    EqualityConstraint operator==(const AffineExpression &lhs, const double &zero);

}




class EcosWrapper {
    size_t n_variables = 0;

    /* Set of named tensor variables in the optimization problem */
    map<string, vector<size_t>> tensor_variable_dimensions;
    map<string, vector<size_t>> tensor_variable_indices;

    /* High-level optimization problem description */
    vector<optimization_problem::SecondOrderConeConstraint> secondOrderConeConstraints;
    vector<optimization_problem::PostiveConstraint> postiveConstraints;
    vector<optimization_problem::EqualityConstraint> equalityConstraints;
    optimization_problem::AffineExpression costFunction;

    /* ECOS problem parameters */
    idxint         ecos_n_variables;
    idxint         ecos_n_constraint_rows;
    idxint         ecos_n_equalities;
    idxint         ecos_n_positive_constraints;
    idxint         ecos_n_cone_constraints;
    vector<idxint> ecos_cone_constraint_dimensions;
    idxint         ecos_n_exponential_cones;
    vector<optimization_problem::Parameter> ecos_G_data_CSS;
    vector<idxint> ecos_G_columns_CCS;
    vector<idxint> ecos_G_rows_CCS;
    vector<optimization_problem::Parameter> ecos_A_data_CCS;
    vector<idxint> ecos_A_columns_CCS;
    vector<idxint> ecos_A_rows_CCS;
    vector<optimization_problem::Parameter> ecos_cost_function_weights;
    vector<optimization_problem::Parameter> ecos_h;
    vector<optimization_problem::Parameter> ecos_b;


    size_t allocate_variable_index();

public:
    void create_tensor_variable(const string &name, const vector<size_t> &dimensions);
    size_t get_tensor_variable_index(const string &name, const vector<size_t> &indices);
    optimization_problem::Variable get_variable(const string &name, const vector<size_t> &indices);

    void add_constraint(optimization_problem::SecondOrderConeConstraint c);
    void add_constraint(optimization_problem::PostiveConstraint c);
    void add_constraint(optimization_problem::EqualityConstraint c);
    void set_cost_function(optimization_problem::AffineExpression c);
    void compile_problem_structure();
};
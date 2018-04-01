#include "EcosWrapper.hpp"
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <tuple>
#include <iostream>

using std::pair;
using std::tuple;
using std::make_pair;








inline size_t tensor_index(const vector<size_t> &indices, const vector<size_t> &dimensions) {
    assert(indices.size() == dimensions.size());
    size_t index = 0;
    for (size_t d = 0; d < indices.size(); ++d) index = index * dimensions[d] + indices[d];
    return index;
}

size_t EcosWrapper::allocate_variable_index() {
    size_t i = n_variables;
    n_variables++;
    return i;
}

void EcosWrapper::create_tensor_variable(const string &name, const vector<size_t> &dimensions) {
    size_t tensor_size = 1;
    for(auto d:dimensions) 
        tensor_size *= d;
    vector<size_t> new_variable_indices(tensor_size);
    for(auto &i:new_variable_indices)
        i = allocate_variable_index();

    tensor_variable_dimensions[name] = dimensions;
    tensor_variable_indices[name] = new_variable_indices;
}

size_t EcosWrapper::get_tensor_variable_index(const string &name, const vector<size_t> &indices) {
    assert(tensor_variable_indices.count(name) > 0);
    auto dims = tensor_variable_dimensions[name];
    assert(indices.size() == dims.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        assert(indices[i] < dims[i]);
    }
    return tensor_variable_indices[name][tensor_index(indices,dims)];
}

optimization_problem::Variable EcosWrapper::get_variable(const string &name, const vector<size_t> &indices) {
    optimization_problem::Variable var;
    var.name = name;
    var.tensor_indices = indices;
    var.problem_index = get_tensor_variable_index(name, indices);
    return var;
}

void EcosWrapper::add_constraint(optimization_problem::SecondOrderConeConstraint c) {
    secondOrderConeConstraints.push_back(c);
}

void EcosWrapper::add_constraint(optimization_problem::PostiveConstraint c) {
    postiveConstraints.push_back(c);
}

void EcosWrapper::add_constraint(optimization_problem::EqualityConstraint c) {
    equalityConstraints.push_back(c);
}

void EcosWrapper::add_minimization_term(optimization_problem::AffineExpression c) {
    costFunction = costFunction + c;
}


template<typename T>
inline bool contains(const vector<T> &v, const T &x)
{ return std::find(v.begin(), v.end(), x) != v.end(); }



bool check_unique_variables_in_affine_expression(const optimization_problem::AffineExpression &affineExpression) {
    // check if a variable is used more than once in an expression

    vector<size_t> variable_indices;
    for(auto const& term:affineExpression.terms) {
        if(term.variable) { // only consider linear terms, not constant terms
            size_t idx = term.variable.value().problem_index;
            if(contains(variable_indices, idx)) {
                // duplicate found!
                return false;
            }
            variable_indices.push_back(idx);
        }
    }
    return true;
}

size_t count_constants_in_affine_expression(const optimization_problem::AffineExpression &affineExpression) {
    size_t count = 0;
    for(auto const& term:affineExpression.terms) {
        if(!term.variable) { // only consider constant terms, not linear terms
            count++;
        }
    }
    return count;
}

void error_check_affine_expression(const optimization_problem::AffineExpression &affineExpression) {
    if(!check_unique_variables_in_affine_expression(affineExpression)) {
        throw std::runtime_error("Error: Duplicate variable in the expression: \n" + affineExpression.print());
    }
    if(count_constants_in_affine_expression(affineExpression) > 1) {
        throw std::runtime_error("Error: More than one constant in the expression: \n" + affineExpression.print());
    }
}

optimization_problem::Parameter get_constant_or_zero(const optimization_problem::AffineExpression &affineExpression) {
    for(auto const& term:affineExpression.terms) {
        if(!term.variable) { // only consider constant terms, not linear terms
            return term.parameter;
        }
    }
    return optimization_problem::Parameter(0.0);
}

// convert sparse matrix format "Dictionary of keys" to "column compressed storage"
void sparse_DOK_to_CCS(
    const map< pair<idxint, idxint>, optimization_problem::Parameter > &sparse_DOK, 
    vector<optimization_problem::Parameter> &data_CCS, 
    vector<idxint> &columns_CCS, 
    vector<idxint> &rows_CCS,
    size_t n_columns) 
{
    using std::get;
    assert(data_CCS.empty());
    assert(columns_CCS.empty());
    assert(rows_CCS.empty());


    // convert to coordinate list
    vector< tuple<idxint, idxint, optimization_problem::Parameter> > sparse_COO;
    for (auto const& e : sparse_DOK) {
        sparse_COO.emplace_back(e.first.first, e.first.second, e.second);
        //               row index ^    column index ^      value ^
    }

    // sort coordinate list by column, then by row
    std::sort(sparse_COO.begin(), sparse_COO.end(), 
        []( const tuple<idxint, idxint, optimization_problem::Parameter> & a, 
            const tuple<idxint, idxint, optimization_problem::Parameter> & b) -> bool
    {
        // define coordinate list order
        if(get<1>(a) == get<1>(b)) {
            return get<0>(a) < get<0>(b);
        }
        return get<1>(a) < get<1>(b);
    });

    // create CCS format
    vector<size_t> non_zero_count_per_column(n_columns, 0);
    for (size_t i = 0; i < sparse_COO.size(); ++i) {
        data_CCS.push_back(get<2>(sparse_COO[i]));
        rows_CCS.push_back(get<0>(sparse_COO[i]));
        non_zero_count_per_column[get<1>(sparse_COO[i])]++;
    }
    columns_CCS.push_back(0);
    for(auto column_count:non_zero_count_per_column) {
        columns_CCS.push_back(column_count + columns_CCS.back());
    }
}

void copy_affine_expression_linear_parts_to_sparse_DOK(
    map< pair<idxint, idxint>, optimization_problem::Parameter > &sparse_DOK, 
    const optimization_problem::AffineExpression &affineExpression,
    size_t row_index) 
{
    for(auto const& term:affineExpression.terms) {
        if(term.variable) { // only consider linear terms, not constant terms
            size_t column_index = term.variable.value().problem_index;
            sparse_DOK[make_pair(row_index, column_index)] = term.parameter;
        }
    }
}

void EcosWrapper::compile_problem_structure() {
    ecos_cone_constraint_dimensions.clear();
    ecos_G_data_CCS.clear();
    ecos_G_columns_CCS.clear();
    ecos_G_rows_CCS.clear();
    ecos_A_data_CCS.clear();
    ecos_A_columns_CCS.clear();
    ecos_A_rows_CCS.clear();
    ecos_cost_function_weights.clear();
    ecos_h.clear();
    ecos_b.clear();


    /* ECOS size parameters */
    ecos_solution_vector.resize(n_variables);
    ecos_n_variables = n_variables;
    ecos_n_cone_constraints = secondOrderConeConstraints.size();
    ecos_n_equalities = equalityConstraints.size();
    ecos_n_positive_constraints = postiveConstraints.size();
    ecos_n_constraint_rows = postiveConstraints.size();
    ecos_n_exponential_cones = 0; // Exponential cones are not supported.
    for(auto const& cone: secondOrderConeConstraints) {
        ecos_n_constraint_rows += 1 + cone.lhs.arguments.size();
        ecos_cone_constraint_dimensions.push_back(1 + cone.lhs.arguments.size());
    }


    /* Error checking for the problem description */
    for(auto const& cone: secondOrderConeConstraints) {
        error_check_affine_expression(cone.rhs);
        for(auto const& affine_expression:cone.lhs.arguments) {
            error_check_affine_expression(affine_expression);
        }
    }
    for(auto const& postiveConstraint:postiveConstraints) {
        error_check_affine_expression(postiveConstraint.lhs);
    }
    for(auto const& equalityConstraint:equalityConstraints) {
        error_check_affine_expression(equalityConstraint.lhs);
    }
    error_check_affine_expression(costFunction);


    /* Build equality constraint parameters (b - A*x == 0) */
    {
        // Construct the sparse A matrix in the "Dictionary of keys" format
        map< pair<idxint, idxint>, optimization_problem::Parameter > A_sparse_DOK;
        vector< optimization_problem::Parameter > b(equalityConstraints.size());

        for (size_t i = 0; i < equalityConstraints.size(); ++i) {
            auto const& affine_expression = equalityConstraints[i].lhs;
            b[i] = get_constant_or_zero(affine_expression);
            copy_affine_expression_linear_parts_to_sparse_DOK(A_sparse_DOK, affine_expression, i);
        }

        // Convert A to "column compressed storage"
        sparse_DOK_to_CCS(A_sparse_DOK, ecos_A_data_CCS, ecos_A_columns_CCS, ecos_A_rows_CCS, ecos_n_variables);
        ecos_b = b;
    }


    /* Build inequality constraint parameters */
    {
        // Construct the sparse G matrix in the "Dictionary of keys" format
        map< pair<idxint, idxint>, optimization_problem::Parameter > G_sparse_DOK;
        vector< optimization_problem::Parameter > h(ecos_n_constraint_rows);

        size_t row_index = 0;

        for(const auto& postiveConstraint:postiveConstraints) {
            h[row_index] = get_constant_or_zero(postiveConstraint.lhs);
            copy_affine_expression_linear_parts_to_sparse_DOK(G_sparse_DOK, postiveConstraint.lhs, row_index);
            row_index++;
        }

        for(const auto& secondOrderConeConstraint:secondOrderConeConstraints) {
            h[row_index] = get_constant_or_zero(secondOrderConeConstraint.rhs);
            copy_affine_expression_linear_parts_to_sparse_DOK(G_sparse_DOK, secondOrderConeConstraint.rhs, row_index);
            row_index++;

            for(const auto& norm2argument:secondOrderConeConstraint.lhs.arguments) {
                h[row_index] = get_constant_or_zero(norm2argument);
                copy_affine_expression_linear_parts_to_sparse_DOK(G_sparse_DOK, norm2argument, row_index);
                row_index++;
            }            
        }

        assert(row_index == ecos_n_constraint_rows); // all rows used?

        // Convert G to "column compressed storage"
        sparse_DOK_to_CCS(G_sparse_DOK, ecos_G_data_CCS, ecos_G_columns_CCS, ecos_G_rows_CCS, ecos_n_variables);
        ecos_h = h;
    }

    /* Build cost function parameters */
    {
        vector< optimization_problem::Parameter > c(ecos_n_variables);
        for(const auto& term:costFunction.terms) {
            if(term.variable) {
                c[term.variable.value().problem_index] = term.parameter;
            }
        }
        ecos_cost_function_weights = c;
    }
}

inline vector<double> get_parameter_values(const vector<optimization_problem::Parameter> &params, double factor) {
    vector<double> result(params.size(), 0.0);
    for (size_t i = 0; i < params.size(); ++i) {
        result[i] = params[i].get_value() * factor;
    }
    return result;
}

void EcosWrapper::solve_problem() {

    vector<double> ecos_cost_function_weights_values = get_parameter_values(ecos_cost_function_weights, 1.0);
    vector<double> ecos_h_values                     = get_parameter_values(ecos_h, 1.0);
    vector<double> ecos_b_values                     = get_parameter_values(ecos_b, 1.0);
    vector<double> ecos_G_data_CCS_values            = get_parameter_values(ecos_G_data_CCS, -1.0);
    vector<double> ecos_A_data_CCS_values            = get_parameter_values(ecos_A_data_CCS, -1.0);
    // The signs for A and G must be flipped because they are negative in the ECOS interface


    pwork *mywork = ECOS_setup(
        ecos_n_variables,
        ecos_n_constraint_rows,
        ecos_n_equalities,
        ecos_n_positive_constraints,
        ecos_n_cone_constraints,
        ecos_cone_constraint_dimensions.data(),
        ecos_n_exponential_cones,
        ecos_G_data_CCS_values.data(),
        ecos_G_columns_CCS.data(),
        ecos_G_rows_CCS.data(),
        ecos_A_data_CCS_values.data(),
        ecos_A_columns_CCS.data(),
        ecos_A_rows_CCS.data(),
        ecos_cost_function_weights_values.data(),
        ecos_h_values.data(),
        ecos_b_values.data()
    );

    if( mywork ){
        ecos_exitflag = ECOS_solve(mywork); 

        // copy solution
        for (int i = 0; i < ecos_n_variables; ++i) {
            ecos_solution_vector[i] = mywork->x[i];
        }
    }
    ECOS_cleanup(mywork, 0); // TODO maybe this does not need to be allocated and freed for every call? Reuse the pwork?
}

void EcosWrapper::print_problem(std::ostream &out) {
    using std::endl;
    out << "Minimize" << endl;
    out << costFunction.print() << endl;

    out << endl << "Subject to equality constraints" << endl;
    for(const auto & equalityConstraint:equalityConstraints) {
        out << equalityConstraint.print() << endl;
    }

    out << endl << "Subject to linear inequalities" << endl;
    for(const auto & postiveConstraint:postiveConstraints) {
        out << postiveConstraint.print() << endl;
    }

    out << endl << "Subject to cone constraints" << endl;
    for(const auto & secondOrderConeConstraint:secondOrderConeConstraints) {
        out << secondOrderConeConstraint.print() << endl;
    }

}
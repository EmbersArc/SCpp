#include "EcosWrapper.hpp"
#include <stdexcept>
#include <algorithm>
#include <map>
#include <utility>
#include <tuple>

using std::map;
using std::pair;
using std::tuple;
using std::make_pair;

namespace optimization_problem {

    string Variable::print() const {
        std::ostringstream s;
        s << name;
        if(tensor_indices.size()>0) {
            s << "[";
            for (size_t i = 0; i < tensor_indices.size(); ++i) {
                if(i) s << ",";
                s << tensor_indices[i];
            }
            s << "]";
        }
        s << "@" << problem_index;
        return s.str();
    }

    string Parameter::print() const {
        std::ostringstream s;
        s << "(" << get_value()  << ")";
        return s.str();
    }

    Parameter::operator AffineTerm() {
        AffineTerm result;
        result.parameter = *this;
        return result;
    }


    Parameter::operator AffineExpression() {
        return (AffineTerm)(*this);
    }


    AffineTerm::operator AffineExpression() { 
        AffineExpression result;
        result.terms.push_back(*this);
        return result;
    }

    string AffineTerm::print() const {
        std::ostringstream s;
        s << parameter.print();
        if(variable) s << "*" << variable.value().print();
        return s.str();
    }

    string AffineExpression::print() const {
        std::ostringstream s;
        for (size_t i = 0; i < terms.size(); ++i) {
            if(i) s << "  + ";
            s << terms[i].print();
        }
        return s.str();
    }

    AffineTerm operator*(const Parameter &parameter, const Variable &variable) {
        AffineTerm affineTerm;
        affineTerm.parameter = parameter;
        affineTerm.variable = variable;
        return affineTerm;
    }


    AffineTerm operator*(const double &const_parameter, const Variable &variable) {
        AffineTerm affineTerm;
        affineTerm.parameter = Parameter(const_parameter);
        affineTerm.variable = variable;
        return affineTerm;        
    }

    AffineExpression operator+(const AffineExpression &lhs, const AffineExpression &rhs) {
        AffineExpression result;
        result.terms.insert(result.terms.end(), lhs.terms.begin(), lhs.terms.end());
        result.terms.insert(result.terms.end(), rhs.terms.begin(), rhs.terms.end());
        return result;
    }

    AffineExpression operator+(const AffineExpression &lhs, const double &rhs) {
        AffineExpression result;
        result.terms.insert(result.terms.end(), lhs.terms.begin(), lhs.terms.end());
        result.terms.push_back(Parameter(rhs));
        return result;
    }


    string Norm2::print() const {
        std::ostringstream s;
        s << "norm2([ ";
        for (size_t i = 0; i < arguments.size(); ++i)
        {
            if(i) s << ", ";
            s << arguments[i].print();
        }
        s << " ])";
        return s.str();
    }

    string PostiveConstraint::print() const {
        std::ostringstream s;
        s << lhs.print() << " >= 0";
        return s.str();
    }

    string EqualityConstraint::print() const {
        std::ostringstream s;
        s << lhs.print() << " == 0";
        return s.str();
    }


    string SecondOrderConeConstraint::print() const {
        std::ostringstream s;
        s << lhs.print() << " <= " << rhs.print();
        return s.str();
    }


    Norm2 norm2(const vector<AffineExpression> &affineExpressions) {
        Norm2 n;
        n.arguments = affineExpressions;
        return n;
    }

    SecondOrderConeConstraint operator<=(const Norm2 &lhs, const AffineExpression &rhs) {
        SecondOrderConeConstraint socc;
        socc.lhs = lhs;
        socc.rhs = rhs;
        return socc;
    }

    SecondOrderConeConstraint operator<=(const Norm2 &lhs, const double &rhs) {
        SecondOrderConeConstraint socc;
        socc.lhs = lhs;
        socc.rhs = Parameter(rhs);
        return socc;
    }

    PostiveConstraint operator>=(const AffineExpression &lhs, const double &zero) {
        assert(zero == 0.0);
        PostiveConstraint result;
        result.lhs = lhs;
        return result;
    }

    EqualityConstraint operator==(const AffineExpression &lhs, const double &zero) {
        assert(zero == 0.0);
        EqualityConstraint result;
        result.lhs = lhs;
        return result;
    }

}










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
    return tensor_variable_indices[name][tensor_index(indices,tensor_variable_dimensions[name])];
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

void EcosWrapper::set_cost_function(optimization_problem::AffineExpression c) {
    costFunction = c;
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

    // create CSS format
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

void EcosWrapper::compile_problem_structure() {
    ecos_cone_constraint_dimensions.clear();
    ecos_G_data_CSS.clear();
    ecos_G_columns_CCS.clear();
    ecos_G_rows_CCS.clear();
    ecos_A_data_CCS.clear();
    ecos_A_columns_CCS.clear();
    ecos_A_rows_CCS.clear();
    ecos_cost_function_weights.clear();
    ecos_h.clear();
    ecos_b.clear();


    /* ECOS size parameters */
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

            for(auto const& term:affine_expression.terms) {
                if(term.variable) { // only consider linear terms, not constant terms
                    size_t column_index = term.variable.value().problem_index;
                    A_sparse_DOK[make_pair(i, column_index)] = term.parameter;
                }
            }
        }

        // Convert A to "column compressed storage"
        sparse_DOK_to_CCS(A_sparse_DOK, ecos_A_data_CCS, ecos_A_columns_CCS, ecos_A_rows_CCS, ecos_n_variables);
        ecos_b = b;
    }
}
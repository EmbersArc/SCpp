#pragma once

#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <sstream> 
#include <experimental/optional> 
using std::vector;
using std::string;
using std::map;
using std::experimental::optional;

namespace optimization_problem {
    struct Variable { // represents an optimization variable x_i
        string name;
        vector<size_t> tensor_indices;
        size_t problem_index;
        string print();
    };

    struct AffineTerm;

    struct AffineExpression { // represents a term like (p_1*x_1 + p_2*x_2 + b)
        vector<AffineTerm> terms;
        string print();
    };

    struct Parameter { // represents a parameter p_i in the opt-problem that can be changed between problem evaluations.
        const double* value_ptr = NULL;
        string print();
        operator AffineTerm();
        operator AffineExpression();
    };

    struct AffineTerm { // represents a linear term (p_1*x_1) or constant term (p_1)
        Parameter parameter;
        optional<Variable> variable;
        string print();
        operator AffineExpression();
    };

    struct Norm2 { // represents a term like norm2([p_1*x_1 + p_2*x_2 + b_1,   p_3*x_3 + p_4*x_4 + b_2 ])
        vector<AffineExpression> arguments;
        string print();
    };

    // represents a term like 
    //      norm2([p_1*x_1 + p_2*x_2 + b_1,   p_3*x_3 + p_4*x_4 + b_2 ])
    //        <= p_5*x_5 + p_6*x_6 + b_3
    struct SecondOrderConeConstraint {
        Norm2 lhs;
        AffineExpression rhs;
        string print();
    };

    AffineTerm operator*(const Parameter &parameter, const Variable &variable);
    AffineExpression operator+(const AffineExpression &lhs, const AffineExpression &rhs);
    Norm2 norm2(const vector<AffineExpression> &affineExpressions);
    SecondOrderConeConstraint operator<=(const Norm2 &lhs, const AffineExpression &rhs);

}


inline size_t tensor_index(const vector<size_t> &indices, const vector<size_t> &dimensions) {
    assert(indices.size() == dimensions.size());
    size_t index = 0;
    for (size_t d = 0; d < indices.size(); ++d) index = index * dimensions[d] + indices[d];
    return index;
}



class EcosWrapper {
    size_t n_variables = 0;

    map<string, vector<size_t>> tensor_variable_dimensions;
    map<string, vector<size_t>> tensor_variable_indices;


    size_t allocate_variable_index() {
        size_t i = n_variables;
        n_variables++;
        return i;
    }

public:
    void create_tensor_variable(const string &name, const vector<size_t> &dimensions) {
        size_t tensor_size = 1;
        for(auto d:dimensions) 
            tensor_size *= d;
        vector<size_t> new_variable_indices(tensor_size);
        for(auto &i:new_variable_indices)
            i = allocate_variable_index();

        tensor_variable_dimensions[name] = dimensions;
        tensor_variable_indices[name] = new_variable_indices;
    }

    size_t get_tensor_variable_index(const string &name, const vector<size_t> &indices) {
        return tensor_variable_indices[name][tensor_index(indices,tensor_variable_dimensions[name])];
    }

    optimization_problem::Variable get_variable(const string &name, const vector<size_t> &indices) {
        optimization_problem::Variable var;
        var.name = name;
        var.tensor_indices = indices;
        var.problem_index = get_tensor_variable_index(name, indices);
        return var;
    }
};
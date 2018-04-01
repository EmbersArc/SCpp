#pragma once

#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <sstream> 
#include <experimental/optional> 
#include <utility>
#include <functional>

using std::vector;
using std::map;
using std::string;
using std::experimental::optional;
using std::pair;


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
        double evaluate(const vector<double> &soln_values) const;
    };

    enum class ParameterSource { Constant, Pointer, Callback };

    // Represents a parameter p_i in the opt-problem that can be changed between problem evaluations.
    // The parameter value can either be constant or dynamically accessed through a pointer or callback.
    class Parameter { 
        std::function<double()> callback;
        const double* dynamic_value_ptr = NULL;
        double const_value = 0;
        ParameterSource parameterSource;
    public:
        explicit Parameter(std::function<double()> callback):callback(callback),parameterSource(ParameterSource::Callback) {
            if(!callback) throw std::runtime_error("Parameter(callback), Invalid Callback Error");
        }
        explicit Parameter(const double* dynamic_value_ptr):dynamic_value_ptr(dynamic_value_ptr),parameterSource(ParameterSource::Pointer) {
            if(dynamic_value_ptr == NULL) throw std::runtime_error("Parameter(NULL), Null Pointer Error");
        }
        explicit Parameter(double const_value):const_value(const_value),parameterSource(ParameterSource::Constant){}
        Parameter():const_value(0),parameterSource(ParameterSource::Constant){}
        double get_value() const {
            if(parameterSource == ParameterSource::Callback) {
                return callback();
            }
            if(parameterSource == ParameterSource::Pointer) {
                return *dynamic_value_ptr;
            }
            return const_value;
        }
        string print() const;
        operator AffineTerm();
        operator AffineExpression();
    };

    struct AffineTerm { // represents a linear term (p_1*x_1) or constant term (p_1)
        Parameter parameter = Parameter(0.0);
        optional<Variable> variable; // a missing Variable represents a constant 1.0
        string print() const;
        double evaluate(const vector<double> &soln_values) const;
        operator AffineExpression();
    };

    struct Norm2 { // represents a term like norm2([p_1*x_1 + p_2*x_2 + b_1,   p_3*x_3 + p_4*x_4 + b_2 ])
        vector<AffineExpression> arguments;
        string print() const;
        double evaluate(const vector<double> &soln_values) const;
    };

    // represents a constraint like 
    //      norm2([p_1*x_1 + p_2*x_2 + b_1,   p_3*x_3 + p_4*x_4 + b_2 ])
    //        <= p_5*x_5 + p_6*x_6 + b_3
    struct SecondOrderConeConstraint {
        Norm2 lhs;
        AffineExpression rhs;
        string print() const;
        double evaluate(const vector<double> &soln_values) const;
    };


    // represents a constraint like 
    //     p_1*x_1 + p_2*x_2 + b >= 0
    struct PostiveConstraint
    {
        AffineExpression lhs;
        string print() const;
        double evaluate(const vector<double> &soln_values) const;
    };

    // represents a constraint like 
    //     p_1*x_1 + p_2*x_2 + b == 0
    struct EqualityConstraint
    {
        AffineExpression lhs;
        string print() const;
        double evaluate(const vector<double> &soln_values) const;
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



    class GenericOptimizationProblem {
    protected:
        size_t n_variables = 0;

        /* Set of named tensor variables in the optimization problem */
        map<string, vector<size_t>> tensor_variable_dimensions;
        map<string, vector<size_t>> tensor_variable_indices;


        size_t allocate_variable_index();
    public:
        void create_tensor_variable(const string &name, const vector<size_t> &dimensions);
        size_t get_tensor_variable_index(const string &name, const vector<size_t> &indices);
        Variable get_variable(const string &name, const vector<size_t> &indices);
        size_t get_n_variables(){ return n_variables; }
    };


    struct SecondOrderConeProgram : public GenericOptimizationProblem {

        vector<SecondOrderConeConstraint> secondOrderConeConstraints;
        vector<PostiveConstraint> postiveConstraints;
        vector<EqualityConstraint> equalityConstraints;
        AffineExpression costFunction = Parameter(0.0);

        void add_constraint(SecondOrderConeConstraint c);
        void add_constraint(PostiveConstraint c);
        void add_constraint(EqualityConstraint c);
        void add_minimization_term(AffineExpression c);
        void print_problem(std::ostream &out);

        bool feasibility_check(const vector<double> &soln_values);
    };
}
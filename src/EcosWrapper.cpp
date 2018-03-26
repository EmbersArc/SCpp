#include "EcosWrapper.hpp"
#include <stdexcept>

namespace optimization_problem {

    string Variable::print() {
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

    string Parameter::print() {
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

    string AffineTerm::print() {
        std::ostringstream s;
        s << parameter.print();
        if(variable) s << "*" << variable.value().print();
        return s.str();
    }

    string AffineExpression::print() {
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


    string Norm2::print() {
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

    string PostiveConstraint::print() {
        std::ostringstream s;
        s << lhs.print() << " >= 0";
        return s.str();
    }

    string EqualityConstraint::print() {
        std::ostringstream s;
        s << lhs.print() << " == 0";
        return s.str();
    }


    string SecondOrderConeConstraint::print() {
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
    // TODO
}

void EcosWrapper::add_constraint(optimization_problem::PostiveConstraint c) {
    // TODO
}

void EcosWrapper::add_constraint(optimization_problem::EqualityConstraint c) {
    // TODO
}

void EcosWrapper::set_cost_function(optimization_problem::AffineExpression c) {
    // TODO
}

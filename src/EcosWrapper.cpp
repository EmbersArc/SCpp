#include "EcosWrapper.hpp"


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
        if(value_ptr == NULL) return "(NULL)";
        std::ostringstream s;
        s << "(" << *value_ptr << ")";
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

    AffineExpression operator+(const AffineExpression &lhs, const AffineExpression &rhs) {
        AffineExpression result;
        result.terms.insert(result.terms.end(), lhs.terms.begin(), lhs.terms.end());
        result.terms.insert(result.terms.end(), rhs.terms.begin(), rhs.terms.end());
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
    };



}
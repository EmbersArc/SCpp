#include <iostream>
#include <string>
#include <memory>
#include "monty.h"
#include "fusion.h"

#include "mosekWrapper.hpp"

using namespace mosek::fusion;
using namespace monty;
using std::make_shared;
using std::shared_ptr;

Expression::t convert_affine_expression(
    const optimization_problem::AffineExpression &expr,
    const vector<Variable::t> &mosek_variables)
{
    shared_ptr<ndarray<Expression::t, 1>> exprs = make_shared<ndarray<Expression::t, 1>>(shape(expr.terms.size()));
    for (size_t i = 0; i < expr.terms.size(); ++i)
    {
        auto &affineTerm = expr.terms[i];

        if (affineTerm.variable)
        {
            (*exprs)[i] = Expr::mul(mosek_variables[affineTerm.variable.value().problem_index], affineTerm.parameter.get_value());
        }
        else
        {
            (*exprs)[i] = Expr::constTerm(1, affineTerm.parameter.get_value());
        }
    }
    return Expr::sum(Expr::vstack(exprs));
}

void MosekWrapper::solve_problem()
{

    Model::t M = new Model();
    auto _M = finally([&]() { M->dispose(); });

    vector<Variable::t> mosek_variables;
    for (size_t i = 0; i < socp.get_n_variables(); ++i)
    {
        mosek_variables.push_back(M->variable());
    }

    for (const auto &equalityConstraint : socp.equalityConstraints)
    {
        M->constraint(convert_affine_expression(equalityConstraint.lhs, mosek_variables), Domain::equalsTo(0.0));
    }

    for (const auto &postiveConstraint : socp.postiveConstraints)
    {
        M->constraint(convert_affine_expression(postiveConstraint.lhs, mosek_variables), Domain::greaterThan(0.0));
    }

    for (const auto &secondOrderConeConstraint : socp.secondOrderConeConstraints)
    {
        const int n = secondOrderConeConstraint.lhs.arguments.size();
        shared_ptr<ndarray<Expression::t, 1>> exprs = make_shared<ndarray<Expression::t, 1>>(shape(n + 1));
        (*exprs)[0] = convert_affine_expression(secondOrderConeConstraint.rhs, mosek_variables);
        for (int i = 0; i < n; ++i)
        {
            (*exprs)[i + 1] = convert_affine_expression(secondOrderConeConstraint.lhs.arguments.at(i), mosek_variables);
        }
        M->constraint(Expr::hstack(exprs), Domain::inQCone(n + 1));
    }

    M->objective(ObjectiveSense::Minimize, convert_affine_expression(socp.costFunction, mosek_variables));

    M->setLogHandler([](const std::string &msg) { print(msg) });
    M->solve();

    solution_vector.resize(socp.get_n_variables());
    for (size_t i = 0; i < socp.get_n_variables(); ++i)
    {
        solution_vector[i] = (*(mosek_variables[i]->level()))[0];
    }
}
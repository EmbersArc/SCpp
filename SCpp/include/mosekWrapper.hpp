#pragma once

#include "optimizationProblem.hpp"

class MosekWrapper
{

    optimization_problem::SecondOrderConeProgram &socp;

    /* result */
    vector<double> solution_vector;

  public:
    explicit MosekWrapper(optimization_problem::SecondOrderConeProgram &_socp) : socp(_socp) {}

    void solveProblem();

    double getSolutionValue(size_t problem_index) const
    {
        return solution_vector[problem_index];
    }

    double getSolutionValue(const string &name, const vector<size_t> &indices)
    {
        return solution_vector[socp.get_variable(name, indices).problem_index];
    }

    vector<double> getSolutionVector()
    {
        if (socp.get_n_variables() > 0 && solution_vector.size() == size_t(socp.get_n_variables()))
        {
            return solution_vector;
        }
        else
        {
            throw std::runtime_error("getSolutionVector(): Solution unavailable.");
        }
    }
};
#pragma once

#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <sstream>
#include <optional>
#include <utility>
#include <functional>

#include "ecos.h"

#include "optimizationProblem.hpp"

using std::map;
using std::optional;
using std::pair;
using std::string;
using std::vector;

void sparse_DOK_to_CCS(
    const map<pair<idxint, idxint>, op::Parameter> &sparse_DOK,
    vector<op::Parameter> &data_CCS,
    vector<idxint> &columns_CCS,
    vector<idxint> &rows_CCS,
    size_t n_columns);

class EcosWrapper
{
    op::SecondOrderConeProgram &socp;

    /* ECOS problem parameters */
    idxint ecos_n_variables;
    idxint ecos_n_constraint_rows;
    idxint ecos_n_equalities;
    idxint ecos_n_positive_constraints;
    idxint ecos_n_cone_constraints;
    vector<idxint> ecos_cone_constraint_dimensions;
    idxint ecos_n_exponential_cones;
    vector<op::Parameter> ecos_G_data_CCS;
    vector<idxint> ecos_G_columns_CCS;
    vector<idxint> ecos_G_rows_CCS;
    vector<op::Parameter> ecos_A_data_CCS;
    vector<idxint> ecos_A_columns_CCS;
    vector<idxint> ecos_A_rows_CCS;
    vector<op::Parameter> ecos_cost_function_weights;
    vector<op::Parameter> ecos_h;
    vector<op::Parameter> ecos_b;

    /* ECOS result */
    vector<double> ecos_solution_vector;
    idxint ecos_exitflag;

  public:
    explicit EcosWrapper(op::SecondOrderConeProgram &_socp);

    void solveProblem(bool verbose = false);

    double getSolutionValue(size_t problem_index) const
    {
        return ecos_solution_vector[problem_index];
    }

    double getSolutionValue(const string &name, const vector<size_t> &indices)
    {
        return ecos_solution_vector[socp.getVariable(name, indices).problem_index];
    }

    vector<double> getSolutionVector()
    {
        if (ecos_n_variables > 0 && ecos_solution_vector.size() == size_t(ecos_n_variables))
        {
            return ecos_solution_vector;
        }
        else
        {
            throw std::runtime_error("getSolutionVector(): Solution unavailable.");
        }
    }
};
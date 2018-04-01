
#include "OptimizationProblem.hpp"
#include <iostream>
using namespace std;


int main() {
    optimization_problem::GenericOptimizationProblem opt_prob;
    opt_prob.create_tensor_variable("X", {14,20});
    opt_prob.create_tensor_variable("U", {3,20});
    opt_prob.create_tensor_variable("Delta", {});

    cout << "Optimization problem index mapping:" << endl << endl;

    for (size_t i = 0; i < 14; ++i)
    {
        for (size_t j = 0; j < 20; ++j)
        {
            size_t index = opt_prob.get_tensor_variable_index("X", {i,j});
            cout << "X[" << i << "," << j << "] is at " << index << endl;
        }
    }


    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 20; ++j)
        {
            size_t index = opt_prob.get_tensor_variable_index("U", {i,j});
            cout << "U[" << i << "," << j << "] is at " << index << endl;
        }
    }


    size_t index = opt_prob.get_tensor_variable_index("Delta", {});
    cout << "Delta is at " << index << endl;
}
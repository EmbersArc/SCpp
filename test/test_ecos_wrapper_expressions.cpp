
#include "EcosWrapper.hpp"
#include <iostream>
using namespace std;
using namespace optimization_problem;


int main() {

    /*
    Toy problem:

    variables [x, y, z]
    parameters [a, g]
    minimize (-x)
    constraints
        (a*x)^2 + y^2 < 1
        y + g > 0
        x + y == z


    Expected solution (if -1 > g):
        infeasible

    Expected solution (if -1 < g < 0):
        x = (1/a) * sqrt(1-g^2)
        y = -g

    Expected solution (if g > 0):
        x = 1/a
        y = 0

    */
    EcosWrapper wrapper;
    wrapper.create_tensor_variable("x", {});
    wrapper.create_tensor_variable("y", {});
    wrapper.create_tensor_variable("z", {});

    auto var = [&](string name){ return wrapper.get_variable(name,{}); };
    auto param = [](double &param_value){ return optimization_problem::Parameter(&param_value); };

    double a = 0.1;
    double g = -0.98;

    wrapper.set_cost_function(  (-1.0) * var("x")  );
    wrapper.add_constraint(  norm2({  param(a) * var("x"), (1.0) * var("y")  }) <= (1.0)  );
    wrapper.add_constraint(  (1.0) * var("x") + (1.0) * var("y") + (-1.0) * var("z") == 0  );
    wrapper.add_constraint(  (1.0) * var("y") + param(g) >= 0  );


}
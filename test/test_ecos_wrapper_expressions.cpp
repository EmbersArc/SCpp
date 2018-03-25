
#include "EcosWrapper.hpp"
#include <iostream>
using namespace std;
using namespace optimization_problem;


int main() {
    EcosWrapper wrapper;
    wrapper.create_tensor_variable("A", {});
    wrapper.create_tensor_variable("B", {});
    wrapper.create_tensor_variable("C", {});
    wrapper.create_tensor_variable("D", {});

    auto var = [&](string name){return wrapper.get_variable(name,{});};
    auto param = [](double &param_value){
        optimization_problem::Parameter p;
        p.value_ptr = &param_value;
        return p;
    };

    double p_A = 2.5;
    double p_B = 0.1;
    double p_C = 5.4;
    double p_D = 3.4;
    double p_k = 1.4;
    double p_h = -0.4;


    auto expr = 
        norm2({
            param(p_A) * var("A") + param(p_B) * var("B") + param(p_k),
            param(p_D) * var("D")
        })
        <= 
        param(p_C) * var("C") + param(p_D) * var("D") + param(p_h);

    cout << expr.print() << endl;

}
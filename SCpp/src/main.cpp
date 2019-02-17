#include "freeFinalTimeAlgorithm.hpp"

int main()
{
    auto model = std::make_shared<Model>();
    freeFinalTimeAlgorithm scSolver(model);
    scSolver.initialize();
    scSolver.solve();
}
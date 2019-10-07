#include "timing.hpp"

using namespace std::chrono;

double tic()
{
    const duration<double, std::milli> s = system_clock::now().time_since_epoch();

    return s.count();
}

double toc(double start)
{
    return tic() - start;
}

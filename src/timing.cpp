#include "timing.hpp"

double tic()
{
    struct timespec t;
    clock_gettime(CLOCK_REALTIME, &t);
    return double(t.tv_sec) * 1000. + double(t.tv_nsec) / 1000000.;
}

double toc(double start)
{
    return tic() - start;
}

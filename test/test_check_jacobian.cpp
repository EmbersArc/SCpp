#include "check_jacobian.h"
#include <iostream>

using std::cout;
using std::endl;

int main() {

    int tries = 0;
    int errors = 0;
    for (int i = 0; i < 1000; ++i)
    {
        tries++;
        if(!check_jacobian()) errors++;
    }

    cout << "tries " << tries << endl;
    cout << "errors " << errors << endl;
    return 0;
}
#include "check_jacobian.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

int main() {

    vector<double> epsilons {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11};
    vector<double> max_absolute_errors(epsilons.size(), 0);
    vector<double> max_relative_errors(epsilons.size(), 0);

    const double random_radius = 1;

    for (int j = 0; j < 100; ++j) {
        for (size_t i = 0; i < epsilons.size(); ++i) {

            double max_absolute_error = 0;
            double max_relative_error = 0;

            check_jacobian( epsilons[i], random_radius, max_absolute_error, max_relative_error );

            max_absolute_errors[i] = fmax(max_absolute_errors[i], max_absolute_error);
            max_relative_errors[i] = fmax(max_relative_errors[i], max_relative_error);
        }
    }


    cout << " |--------------------|--------------------|--------------------| " << endl;
    cout << " |       epsilon      |   absolute error   |   relative error   | " << endl;
    cout << " |--------------------|--------------------|--------------------| " << endl;
    for (size_t i = 0; i < epsilons.size(); ++i) {
        cout << " | " << setfill(' ') << setw(18) << epsilons[i];
        cout << " | " << setfill(' ') << setw(18) << max_absolute_errors[i];
        cout << " | " << setfill(' ') << setw(18) << max_relative_errors[i];
        cout << " | " << endl;
    }
    cout << " |--------------------|--------------------|--------------------| " << endl;


    return 0;
}
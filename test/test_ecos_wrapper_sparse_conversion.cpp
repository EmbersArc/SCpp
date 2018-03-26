
#include "EcosWrapper.hpp"
#include <iostream>
using namespace std;


int main() {


    map< pair<idxint, idxint>, optimization_problem::Parameter > sparse_DOK;

    /* 
    Test example

        0 0 0 4
        1 0 3 5
        0 0 2 0
        6 0 0 9

    Expected output

    data = 1 6 3 2 4 5 9
    rows = 1 3 1 2 0 1 3
    cols = 0 2 2 4 7
        */

    sparse_DOK[make_pair(1,0)] = optimization_problem::Parameter(1.0);
    sparse_DOK[make_pair(2,2)] = optimization_problem::Parameter(2.0);
    sparse_DOK[make_pair(1,2)] = optimization_problem::Parameter(3.0);
    sparse_DOK[make_pair(0,3)] = optimization_problem::Parameter(4.0);
    sparse_DOK[make_pair(1,3)] = optimization_problem::Parameter(5.0);
    sparse_DOK[make_pair(3,0)] = optimization_problem::Parameter(6.0);
    sparse_DOK[make_pair(3,3)] = optimization_problem::Parameter(9.0);



    vector<optimization_problem::Parameter> data_CCS;
    vector<idxint> columns_CCS;
    vector<idxint> rows_CCS;

    sparse_DOK_to_CCS(sparse_DOK, data_CCS, columns_CCS, rows_CCS, 4);

    cout << "data  ";
    for(auto d:data_CCS) 
        cout << d.get_value() << " ";
    cout << endl;

    cout << "rows  ";
    for(auto d:rows_CCS) 
        cout << d << " ";
    cout << endl;

    cout << "cols  ";
    for(auto d:columns_CCS) 
        cout << d << " ";
    cout << endl;

}
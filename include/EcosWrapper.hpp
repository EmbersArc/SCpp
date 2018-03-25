#pragma once

#include <vector>
#include <map>
#include <string>
#include <cassert>
using std::vector;
using std::string;
using std::map;



inline size_t tensor_index(const vector<size_t> &indices, const vector<size_t> &dimensions) {
    assert(indices.size() == dimensions.size());
    size_t index = 0;
    for (size_t d = 0; d < indices.size(); ++d) index = index * dimensions[d] + indices[d];
    return index;
}

class EcosWrapper {
    size_t n_variables = 0;

    map<string, vector<size_t>> tensor_variable_dimensions;
    map<string, vector<size_t>> tensor_variable_indices;


    size_t allocate_variable_index() {
        size_t i = n_variables;
        n_variables++;
        return i;
    }

public:
    void create_tensor_variable(const string &name, const vector<size_t> &dimensions) {
        size_t tensor_size = 1;
        for(auto d:dimensions) 
            tensor_size *= d;
        vector<size_t> new_variable_indices(tensor_size);
        for(auto &i:new_variable_indices)
            i = allocate_variable_index();

        tensor_variable_dimensions[name] = dimensions;
        tensor_variable_indices[name] = new_variable_indices;
    }

    size_t get_tensor_variable_index(const string &name, const vector<size_t> &indices) {
        return tensor_variable_indices[name][tensor_index(indices,tensor_variable_dimensions[name])];
    }
};
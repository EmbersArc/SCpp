#include <vector>

template <typename T>
std::vector<T> reduce_vector(const std::vector<T> &v, size_t steps)
{
    const size_t size = v.size();

    std::vector<T> new_vector;

    for (size_t i = 0; i < steps; i++)
    {
        const size_t index = size_t(size / steps * i);
        new_vector.push_back(v.at(index));
    }
    return new_vector;
}
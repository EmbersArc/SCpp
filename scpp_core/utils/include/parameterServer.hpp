
#pragma once

#include <string>

#include <Eigen/Dense>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

/**
 * An auxiliary function which loads an Eigen matrix or a scalar from a file.
 * The file uses property tree data structure with INFO format.
 *
 * It has the following format:	<br>
 * 
 * scalarName   value   <br>
 *
 * matrixName	<br>
 * {	<br>
 *   scaling 1e+0				<br>
 *   (0,0) value    ; M(0,0)	<br>
 *   (1,0) value    ; M(1,0)	<br>
 *   (0,1) value    ; M(0,1)	<br>
 *   (1,1) value    ; M(1,1)	<br>
 * } 	<br>
 *
 * If a value for a specific element is not defined it will set by default to zero.
 *
 * For vectors the second index can be left out.
 */
class ParameterServer
{
public:
    explicit ParameterServer(const std::string &filename)
    {
        try
        {
            boost::property_tree::read_info(filename, pt);
        }
        catch (const boost::property_tree::info_parser_error &e)
        {
            fmt::print("Could not open file for reading: {}\n", filename);
            fmt::print("{}\n", e.what());
        }
    }

    template <typename T>
    void loadScalar(
        const std::string &scalarName,
        T &scalar);

    template <typename T>
    void loadMatrix(
        const std::string &matrixName,
        Eigen::MatrixBase<T> &matrix);

private:
    boost::property_tree::ptree pt;
};

template <typename T>
void ParameterServer::loadScalar(
    const std::string &scalarName,
    T &scalar)
{
    try
    {
        scalar = pt.get<T>(scalarName);
    }
    catch (...)
    {
        throw std::runtime_error(fmt::format("WARNING: Failed to load scalar type: {}!\n", scalarName));
    }
}

template <typename T>
void ParameterServer::loadMatrix(
    const std::string &matrixName,
    Eigen::MatrixBase<T> &matrix)
{
    using scalar_t = typename Eigen::MatrixBase<T>::Scalar;

    const scalar_t scaling = pt.get<scalar_t>(matrixName + ".scaling", 1);

    matrix.setZero();

    const size_t rows = matrix.rows();
    const size_t cols = matrix.cols();

    boost::property_tree::ptree matrix_pt = pt.get_child(matrixName);

    const size_t num_entries = matrix_pt.size() - matrix_pt.count("scaling");
    if (num_entries < size_t(matrix.size()))
    {
        throw std::runtime_error(fmt::format("Missing entries in matrix type: {}!\n", matrixName));
    }
    if (num_entries > size_t(matrix.size()))
    {
        throw std::runtime_error(fmt::format("Redundant entries in matrix type: {}!\n", matrixName));
    }

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            try
            {
                if (cols == 1)
                {
                    matrix(i) = matrix_pt.get<scalar_t>(fmt::format("({})", i));
                }
                else
                {
                    matrix(i, j) = matrix_pt.get<scalar_t>(fmt::format("({},{})", i, j));
                }
            }
            catch (...)
            {
                throw std::runtime_error(fmt::format("Failed to load matrix type: {}!\n", matrixName));
            }
        }
    }

    matrix *= scaling;
}
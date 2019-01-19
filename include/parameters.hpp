
#pragma once

#include <string>

#include <Eigen/Dense>
#include <fmt/format.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

using fmt::format;
using fmt::print;
using std::string;

template <typename T>
inline void loadScalar(
    const string &filename,
    const string &scalarName,
    T &scalar)
{
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);

    try
    {
        scalar = pt.get<T>(scalarName);
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error(format("WARNING: Failed to load scalar type: {}!\n", scalarName));
    }
}

/**
 * An auxiliary function which loads an Eigen matrix from a file. The file uses property tree data structure with INFO format (refer to www.goo.gl/fV3yWA).
 *
 * It has the following format:	<br>
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
 * @param [in] filename: File name which contains the configuration data.
 * @param [in] matrixName: The key name assigned to the matrix in the config file.
 * @param [out] matrix: The loaded matrix.
 */
template <typename Derived>
inline void loadMatrix(
    const string &filename,
    const string &matrixName,
    Eigen::MatrixBase<Derived> &matrix)
{
    typedef typename Eigen::MatrixBase<Derived>::Scalar scalar_t;

    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);

    scalar_t scaling = pt.get<scalar_t>(matrixName + ".scaling", 1);

    matrix.setZero();

    size_t rows = matrix.rows();
    size_t cols = matrix.cols();

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            try
            {
                matrix(i, j) = pt.get<scalar_t>(format("{}.({},{})", matrixName, i, j));
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(format("WARNING: Failed to load matrix type: {}!\n", matrixName));
            }
        }
    }

    matrix *= scaling;
}

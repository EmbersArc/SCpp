
#pragma once

#include <string>

#include <Eigen/Dense>
#include <fmt/format.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

using fmt::format;
using fmt::print;
using std::string;

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
    explicit ParameterServer(const string &filename);

    template <typename T>
    void loadScalar(
        const string &scalarName,
        T &scalar);

    template <typename T>
    void loadMatrix(
        const string &matrixName,
        Eigen::MatrixBase<T> &matrix);

  private:
    boost::property_tree::ptree pt;
};

template <typename T>
void ParameterServer::loadScalar(
    const string &scalarName,
    T &scalar)
{
    try
    {
        scalar = pt.get<T>(scalarName);
    }
    catch (...)
    {
        throw std::runtime_error(format("WARNING: Failed to load scalar type: {}!\n", scalarName));
    }
}

template <typename T>
void ParameterServer::loadMatrix(
    const string &matrixName,
    Eigen::MatrixBase<T> &matrix)
{
    typedef typename Eigen::MatrixBase<T>::Scalar scalar_t;

    const scalar_t scaling = pt.get<scalar_t>(matrixName + ".scaling", 1);

    matrix.setZero();

    const size_t rows = matrix.rows();
    const size_t cols = matrix.cols();

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            try
            {
                if (cols == 1)
                {
                    matrix(i) = pt.get<scalar_t>(format("{}.({})", matrixName, i));
                }
                else
                {
                    matrix(i, j) = pt.get<scalar_t>(format("{}.({},{})", matrixName, i, j));
                }
            }
            catch (...)
            {
                throw std::runtime_error(format("WARNING: Failed to load matrix type: {}!\n", matrixName));
            }
        }
    }

    matrix *= scaling;
}
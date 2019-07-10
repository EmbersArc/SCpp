
#include "parameterServer.hpp"

ParameterServer::ParameterServer(const std::string &filename)
{
    try
    {
        boost::property_tree::read_info(filename, pt);
    }
    catch (...)
    {
        fmt::print("Could not open file for reading: {}", filename);
    }
}

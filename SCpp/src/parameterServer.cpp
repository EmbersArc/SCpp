
#include "parameterServer.hpp"

ParameterServer::ParameterServer(const std::string &filename)
{
    boost::property_tree::read_info(filename, pt);
}

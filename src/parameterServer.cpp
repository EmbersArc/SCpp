
#include "parameterServer.hpp"

ParameterServer::ParameterServer(const string &filename)
{
    boost::property_tree::read_info(filename, pt);
}

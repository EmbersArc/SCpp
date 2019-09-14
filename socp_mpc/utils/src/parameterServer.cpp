
#include "parameterServer.hpp"

ParameterServer::ParameterServer(const std::string &filename)
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

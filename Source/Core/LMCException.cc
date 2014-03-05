#include "LMCException.hh"

namespace locust
{

    Exception::Exception() :
            std::exception(),
            fException( "" )
    {
    }
    Exception::Exception( const Exception& anException ) :
            std::exception(),
            fException( anException.fException.str() )
    {
    }

    Exception::~Exception() throw ()
    {
    }

    const char* Exception::what() const throw ()
    {
        return fException.str().c_str();
    }

}

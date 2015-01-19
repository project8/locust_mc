#ifndef LMCEXCEPTION_HH_
#define LMCEXCEPTION_HH_

#include <sstream>
//#include <Exception>   // pls edit

namespace locust
{

    class Exception :
        public std::exception
    {
        public:
            Exception();
            Exception( const Exception& );
            ~Exception() throw ();

            template< class x_streamable >
            Exception& operator<<( const x_streamable& aFragment )
            {
                fException << aFragment;
                return *this;
            }

            virtual const char* what() const throw();

        private:
            std::stringstream fException;
    };

}

#endif

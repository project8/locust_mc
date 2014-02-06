#ifndef LMCSINGLETON_H_
#define LMCSINGLETON_H_

#include "LMCDestroyer.hh"

#include <cstddef>

namespace locust
{

    template< class XType >
    class Singleton
    {
        public:
            static XType* GetInstance();

        private:
            static XType* fInstance;
            static Destroyer< XType > fDestroyer;

        protected:
            Singleton();

            friend class Destroyer< XType >;
            ~Singleton();
    };

    template< class XType >
    XType* Singleton< XType >::fInstance = NULL;

    template< class XType >
    Destroyer< XType > Singleton< XType >::fDestroyer;

    template< class XType >
    XType* Singleton< XType >::GetInstance()
    {
        if( fInstance == NULL )
        {
            fInstance = new XType();
            fDestroyer.SetDoomed( fInstance );
        }
        return fInstance;
    }

    template< class XType >
    Singleton< XType >::Singleton()
    {
    }
    template< class XType >
    Singleton< XType >::~Singleton()
    {
    }

} // namespace locust

#endif /* LMCSINGLETON_H_ */

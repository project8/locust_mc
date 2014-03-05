/*
 * LMCDestroyer.hh
 *
 *  Created on: Nov 7, 2011
 *      Author: nsoblath
 */

#ifndef LMCDESTROYER_H_
#define LMCDESTROYER_H_

namespace locust
{

    template< class XDoomed >
    class Destroyer
    {
        public:
            Destroyer( XDoomed* = 0 );
            ~Destroyer();

            void SetDoomed( XDoomed* );

        private:
            // Prevent users from making copies of a Destroyer to avoid double deletion:
            Destroyer( const Destroyer< XDoomed >& );
            void operator=( const Destroyer< XDoomed >& );

        private:
            XDoomed* fDoomed;
    };

    template< class XDoomed >
    Destroyer< XDoomed >::Destroyer( XDoomed* d )
    {
        fDoomed = d;
    }

    template< class XDoomed >
    Destroyer< XDoomed >::~Destroyer()
    {
        delete fDoomed;
    }

    template< class XDoomed >
    void Destroyer< XDoomed >::SetDoomed( XDoomed* d )
    {
        fDoomed = d;
    }

} /* namespace locust */
#endif /* LMCDESTROYER_H_ */

/*
 * LMCRootTreeWriter.hh
 *
 *  Created on: Mar. 21, 2020
 *      Author: pslocum
 */

#ifndef LMCROOTTREEWRITER_HH_
#define LMCROOTTREEWRITER_HH_

#include "LMCFileWriter.hh"
#include <string>


namespace locust
{
 /*!
 @class RootTreeWriter
 @author P. Slocum
 @brief Derived class to configure the Root TTree that will be written to file.
 @details
 Available configuration options:
 No input parameters
 */
    class RootTreeWriter : public FileWriter, public scarab::singleton< locust::RootTreeWriter >
    {

	public:

    	RootTreeWriter();
        virtual ~RootTreeWriter();

        virtual bool Configure( const scarab::param_node& aNode );
        allow_singleton_access( RootTreeWriter );

        virtual void WriteEvent(Event* anEvent);
        virtual void WriteRunParameters( RunParameters* aRunParameter, const char* aParameterName );

        private:


    };


} /* namespace locust */

#endif /* LMCROOTTREEWRITER_HH_ */

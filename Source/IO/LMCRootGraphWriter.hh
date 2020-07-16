/*
 * LMCRootGraphWriter.hh
 *
 *  Created on: Jul. 10, 2020
 *      Author: pslocum
 */

#ifndef LMCROOTGRAPHWRITER_HH_
#define LMCROOTGRAPHWRITER_HH_

#include "LMCFileWriter.hh"
#include <string>
#include <vector>


namespace locust
{
 /*!
 @class RootGraphWriter
 @author P. Slocum
 @brief Derived class to configure the Root TGraph that will be written to file.
 @details
 Available configuration options:
 No input parameters
 */
    class RootGraphWriter : public FileWriter, public scarab::singleton< locust::RootGraphWriter >
    {

	public:

    	RootGraphWriter();
        virtual ~RootGraphWriter();

        virtual bool Configure( const scarab::param_node& aNode );
        allow_singleton_access( RootGraphWriter );

    private:

        void Write2DGraph(TGraph* aGraph);
        void WriteVectorGraph(std::vector<double> xVector, std::vector<double> yVector);


    };


} /* namespace locust */

#endif /* LMCROOTGRAPHWRITER_HH_ */

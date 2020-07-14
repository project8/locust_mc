/*
 * LMCRootHistoWriter.hh
 *
 *  Created on: Jul. 8, 2020
 *      Author: pslocum
 */

#ifndef LMCROOTHISTOWRITER_HH_
#define LMCROOTHISTOWRITER_HH_

#include "LMCFileWriter.hh"
#include <string>
#include <vector>


namespace locust
{
 /*!
 @class RootHistoWriter
 @author P. Slocum
 @brief Derived class to configure the Root histogram that will be written to file.
 @details
 Available configuration options:
 No input parameters
 */
    class RootHistoWriter : public FileWriter, public scarab::singleton< locust::RootHistoWriter >
    {

	public:

    	RootHistoWriter();
        virtual ~RootHistoWriter();

        virtual bool Configure( const scarab::param_node& aNode );
        allow_singleton_access( RootHistoWriter );

        private:

        virtual void Write1DHisto(TH1D* aHisto);
        virtual void Write2DHisto(TH2D* aHisto);
        virtual void WriteVector1DHisto(std::vector<double> aVector, double xmin, double xmax);


    };


} /* namespace locust */

#endif /* LMCROOTHISTOWRITER_HH_ */

/*
 * LMCRootTreeWriter.hh
 *
 *  Created on: Mar. 21, 2020
 *      Author: pslocum
 */

#ifndef LMCROOTTREEWRITER_HH_
#define LMCROOTTREEWRITER_HH_

#include "TFile.h"  // order of includes matters.
#include "TTree.h"  // include these first.

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

        virtual void WriteRootFile(Event* anEvent);
        virtual void WriteRunParameters( RunParameters* aRunParameter, const char* aParameterName );

        virtual double GetTestVar();
        virtual void SetTestVar(double aValue);

        virtual void OpenFile(std::string aFileFlag);
        virtual void CloseFile();

        private:

        double fTestVar;
        TFile* fFile;
        std::string fRoot_filename;


    };


} /* namespace locust */

#endif /* LMCROOTTREEWRITER_HH_ */

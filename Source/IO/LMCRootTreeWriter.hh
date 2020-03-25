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
#include "LMCEvent.hh"

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
//    class RootTreeWriter: public FileWriter
    class RootTreeWriter : FileWriter
    {
	public:

        RootTreeWriter();
        virtual ~RootTreeWriter();

        virtual bool Configure( const scarab::param_node& aNode );

        void WriteRootFile(Event* anEvent);

        double GetTestVar();

        void SetTestVar(double aValue);

        void OpenFile(std::string aFileName);
        void CloseFile();

        private:
//        double fTestVar;
        TFile* fFile;


    };


} /* namespace locust */

#endif /* LMCROOTTREEWRITER_HH_ */

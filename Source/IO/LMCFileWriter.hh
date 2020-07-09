/*
 * LMCFileWriter.hh
 *
 *  Created on: Mar 21, 2020
 *      Author: pslocum
 */

#ifndef LMCFILEWRITER_HH_
#define LMCFILEWRITER_HH_

#include "TFile.h"  // order of includes matters.
#include "TTree.h"  // include these first.
#include "LMCEvent.hh"
#include "LMCRunParameters.hh"

#include "param.hh"
#include "logger.hh"
#include "singleton.hh"

namespace locust
{
 /*!
 @class FileWriter
 @author P. Slocum
 @brief Base class to contain LMCFileWriter subclasses, such as LMCRootTreeFileWriter.
 @details
 Available configuration options:
 No input parameters
 */

//    class FileWriter : public scarab::singleton< locust::FileWriter >
    class FileWriter
    {
    	protected:
    	FileWriter();
    	virtual ~FileWriter();
        virtual bool Configure( const scarab::param_node& aNode );

    	public:
        virtual double GetTestVar() {};
        virtual void SetTestVar(double aValue) {    	printf("filewriter says hello\n");
};
        virtual void WriteRootFile(Event* anEvent) {};
        virtual void WriteRunParameters( RunParameters* aRunParameter, const char* aParameterName ) {};
        virtual void OpenFile(std::string aFileFlag) {};
        virtual void CloseFile() {};


    	private:

//        double fTestVar;


    };


} /* namespace locust */

#endif /* LMCFILEWRITER_HH_ */

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
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
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
        virtual void SetTestVar(double aValue) {};
        virtual void WriteEvent(Event* anEvent) {};
        virtual void WriteRunParameters( RunParameters* aRunParameter, const char* aParameterName ) {};
        virtual void Write1DHisto(TH1D* aHisto) {};
        virtual void Write2DHisto(TH2D* aHisto) {};
        virtual void WriteVector1DHisto(std::vector<double> aVector, double xmin, double xmax) {};
        virtual void Write2DGraph(TGraph* aGraph) {};
        virtual void WriteVectorGraph(std::vector<double> xVector, std::vector<double> yVector) {};

        virtual void OpenFile(std::string aFileFlag);
        virtual void CloseFile();

        virtual void SetFilename(std::string aFilename);
        virtual std::string GetFilename();


    	private:
        TFile* fFile;
        std::string fRoot_filename;


    };


} /* namespace locust */

#endif /* LMCFILEWRITER_HH_ */

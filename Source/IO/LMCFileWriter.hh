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

    class FileWriter : scarab::singleton< locust::FileWriter >
    {
    	protected:
    	FileWriter();
    	virtual ~FileWriter();
        allow_singleton_access( FileWriter );


    };


} /* namespace locust */

#endif /* LMCFILEWRITER_HH_ */

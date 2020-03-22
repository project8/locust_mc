/*
 * LMCFileWriter.hh
 *
 *  Created on: Mar 21, 2020
 *      Author: pslocum
 */

#ifndef LMCFILEWRITER_HH_
#define LMCFILEWRITER_HH_

#include "param.hh"
#include "logger.hh"

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
    class FileWriter
    {

        public:
            FileWriter();
            virtual ~FileWriter();
            virtual bool Configure( const scarab::param_node& );

        private:
    };


} /* namespace locust */

#endif /* LMCFILEWRITER_HH_ */

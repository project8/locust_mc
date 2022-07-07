/*
 * LMCAnalyticResponseFunction.hh
 *
 *  Created on: Jul 6, 2022
 *      Author: pslocum
 */

#ifndef LMCANALYTICRESPONSEFUNCTION_HH_
#define LMCANALYTICRESPONSEFUNCTION_HH_
#include "param.hh"


namespace locust
{
 /*!
 @class LMCAnalyticResponseFunction
 @author P. Slocum
 @brief Base class to define analytic response functions
 @details
 Available configuration options:
 No input parameters
 */
    class AnalyticResponseFunction
    {

        public:
            AnalyticResponseFunction();
            virtual ~AnalyticResponseFunction();

            virtual bool Configure( const scarab::param_node& aNode );


};


} /* namespace locust */

#endif

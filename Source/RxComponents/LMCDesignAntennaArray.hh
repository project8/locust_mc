/*
 * LMCSDesignAntennaArray.hh
 *
 *  Created on: Oct 29, 2020
 *      Author: Arina Telles
 */

#ifndef LMCDESIGNANTENNAARRAY_HH_
#define LMCDESIGNANTENNAARRAY_HH_

#include "LMCPowerCombiner.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class DesignAntennaArray
 @author A. Telles
 @brief Derived class describing the array of the design antenna
 @details
 Available configuration options:
 No input parameters
 */


    class DesignAntennaArray: public PowerCombiner
    {

        public:
            DesignAntennaArray();
            virtual ~DesignAntennaArray();
            virtual bool Configure( const scarab::param_node& aNode );
            virtual bool SetVoltageDampingFactors();

            virtual Receiver* ChooseElement();


        private:
        	double fImpedanceTransformation;  // != 1.0 either in Locust or in HFSS, but not both.



    };


} /* namespace locust */

#endif /* LMCDESIGNANTENNAARRAY_HH_ */

/*
 * LMCUnitCell.hh
 *
 *  Created on: Feb 28, 2020
 *      Author: pslocum
 */

#ifndef LMCUNITCELL_HH_
#define LMCUNITCELL_HH_

#include "LMCPowerCombiner.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class UnitCell
 @author P. Slocum
 @brief Derived class describing a unit cell.
 @details
 Available configuration options:
 No input parameters
 */


    class UnitCell: public PowerCombiner
    {

        public:
            UnitCell();
            virtual ~UnitCell();
            virtual bool Configure( const scarab::param_node& aNode );
        	virtual bool SetVoltageDampingFactors();



        private:





    };


} /* namespace locust */

#endif /* LMCUNITCELL_HH_ */

/*
 * LMCAliasingUtility.hh
 *
 *  Created on: Feb 10, 2023
 *      Author: pslocum
 */

#ifndef LMCALIASINGUTILITY_HH_
#define LMCALIASINGUTILITY_HH_

#include <math.h>
#include "LMCUtility.hh"

namespace locust
{

    /*!
     @class AliasingUtility
     @author P. L. Slocum Feb. 11 2023

     @brief Class to calculated aliased HF frequencies.

     @details
     Operates in time space

     Configuration name: N/A

     Available configuration options: N/A

     */
    class AliasingUtility : public Utility
    {
    public:

        AliasingUtility();
        virtual ~AliasingUtility();

        bool Configure();
        bool CheckAliasing( double RF, double LO, double fs, double dr );


    private:

    };


} /* namespace locust */

#endif /* LMCALIASINGUTILITY_HH_ */

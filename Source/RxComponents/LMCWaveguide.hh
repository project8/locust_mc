/*
 * LMCWaveguide.hh
 *
 *  Created on: Nov 15, 2021
 *      Author: atelles
 */

#ifndef LMCWAVEGUIDE_HH_
#define LMCWAVEGUIDE_HH_

#include "LMCPowerCombiner.hh"
#include "param.hh"
#include "logger.hh"
#include <iostream>

namespace locust
{
 /*!
 @class Waveguide
 @author A. B. Telles
 @brief Derived class describing the waveguide e-gun apparatus
 @details
 Available configuration options:
 - "impedance-transformation": float -- Transform from transfer function impedance to 50 ohms, sqrt(50/Z_tf)
 */


    class Waveguide: public PowerCombiner
    {

        public:
            Waveguide();
            virtual ~Waveguide();
            virtual bool Configure( const scarab::param_node& aNode );

        private:
        	double fImpedanceTransformation; 



    };


} /* namespace locust */

#endif /* LMCWAVEGUIDE_HH_ */

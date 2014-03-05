/*
 * Constants.hh
 *
 *  Created on: Feb 3, 2014
 *      Author: nsoblath
 */

#ifndef CONSTANTS_HH_
#define CONSTANTS_HH_

namespace locust
{

    struct Constants
    {
        double kB=1.3806e-23;  // Boltzman constant W/K
        double c=3e10; //speed of light in cm/s
        double ecyclo=1.758820e11; //electron cyclotron frequency in rad/s
        double emass=510998.9; //electron mass in eV
        double JoulesToEv=1/1.6e-19; // eV/J
    };

} /* namespace locust */

#endif /* CONSTANTS_HH_ */

/*
 * LMCPozarCylindricalCavity.hh
 *
 *  Created on: Apr. 11, 2023
 *      Author: pslocum
 */

#ifndef LMCPOZARCYLINDRICALCAVITY_HH_
#define LMCPOZARCYLINDRICALCAVITY_HH_

#include "param.hh"

#include "logger.hh"
#include "LMCField.hh"
#include <boost/math/special_functions/bessel.hpp>

#include <vector>

namespace locust
{
 /*!
 @class PozarCylindricalCavity
 @author P. Slocum
 @brief Derived class to define CylindricalCavity fields as in Pozar.
 @details
 Available configuration options:
 No input parameters
 */


    class PozarCylindricalCavity: public FieldCore
    {
        public:

    	    PozarCylindricalCavity();
		    virtual ~PozarCylindricalCavity();

            virtual std::vector<double> TE_E(double R, double L, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TE_H(double R, double L, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TM_E(double R, double L, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TM_H(double R, double L, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);

    };


}; /* namespace locust */

#endif /* LMCPOZARCYLINDRICALCAVITY_HH_ */

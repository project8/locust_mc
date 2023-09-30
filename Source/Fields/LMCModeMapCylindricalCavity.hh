/*
 * LMCModeMapCylindricalCavity.hh
 *
 *  Created on: Sept. 29, 2023
 *      Author: pslocum
 */

#ifndef LMCMODEMAPCYLINDRICALCAVITY_HH_
#define LMCMODEMAPCYLINDRICALCAVITY_HH_

#include "param.hh"

#include "logger.hh"
#include "LMCField.hh"

#include <vector>

namespace locust
{
 /*!
 @class ModeMapCylindricalCavity
 @author P. Slocum
 @brief Derived class to define CylindricalCavity fields from an uploaded mode map file.
 @details
 Available configuration options:
 No input parameters
 */


    class ModeMapCylindricalCavity: public FieldCore
    {
        public:

    	    ModeMapCylindricalCavity();
		    virtual ~ModeMapCylindricalCavity();

            virtual std::vector<double> TE_E(double R, double L, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TE_H(double R, double L, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TM_E(double R, double L, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TM_H(double R, double L, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);

            virtual bool ReadModeMapTE_E(std::string aFilename);

        private:

            std::vector<std::vector<std::vector<double>>> fModeMapTE_E;

    };


}; /* namespace locust */

#endif /* LMCMODEMAPCYLINDRICALCAVITY_HH_ */

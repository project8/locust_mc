/*
 * LMCModeMapCavity.hh
 *
 *  Created on: Sept. 29, 2023
 *      Author: pslocum
 */

#ifndef LMCMODEMAPCAVITY_HH_
#define LMCMODEMAPCAVITY_HH_

#include "param.hh"

#include "logger.hh"
#include "LMCField.hh"

#include <vector>

namespace locust
{
 /*!
 @class ModeMapCavity
 @author P. Slocum
 @brief Derived class to define Cavity fields from an uploaded mode map file.
 @details
 Available configuration options:
 No input parameters
 */


    class ModeMapCavity: public FieldCore
    {
        public:

    	    ModeMapCavity();
		    virtual ~ModeMapCavity();

            virtual std::vector<double> TE_E(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TE_H(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TM_E(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);
            virtual std::vector<double> TM_H(double dim1, double dim2, double dim3, int l, int m, int n, double r, double theta, double z, bool includeOtherPols);

            virtual bool ReadModeMapTE_E(std::string aFilename);

        private:

            std::vector<std::vector<double>> fModeMapTE_E;

    };


}; /* namespace locust */

#endif /* LMCMODEMAPCAVITY_HH_ */

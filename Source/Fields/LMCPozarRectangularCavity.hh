/*
 * LMCPozarRectangularCavity.hh
 *
 *  Created on: Nov. 16, 2023
 *      Author: pslocum
 */

#ifndef LMCPOZARRECTANGULARCAVITY_HH_
#define LMCPOZARRECTANGULARCAVITY_HH_

#include "LMCField.hh"

#include <vector>

namespace locust
{
 /*!
 @class PozarRectangularCavity
 @author P. Slocum
 @brief Derived class to define RectangularCavity fields as in Pozar.
 @details
 Available configuration options:

 */

    class PozarRectangularCavity: public FieldCore
    {
        public:

	        PozarRectangularCavity();
	        virtual ~PozarRectangularCavity();

            virtual std::vector<double> TE_E(double dimX, double dimY, double dimZ, int l, int m, int n, double xKass, double yKass, double zKass, bool includeOtherPols);
            virtual std::vector<double> TE_H(double dimX, double dimY, double dimZ, int l, int m, int n, double xKass, double yKass, double zKass, bool includeOtherPols);
            virtual std::vector<double> TM_E(double dimX, double dimY, double dimZ, int l, int m, int n, double xKass, double yKass, double zKass, bool includeOtherPols);
            virtual std::vector<double> TM_H(double dimX, double dimY, double dimZ, int l, int m, int n, double xKass, double yKass, double zKass, bool includeOtherPols);

    };


}; /* namespace locust */

#endif /* LMCPOZARRECTANGULARCAVITY_HH_ */

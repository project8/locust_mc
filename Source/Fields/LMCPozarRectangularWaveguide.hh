/*
 * LMCPozarRectangularWaveguide.hh
 *
 *  Created on: Apr. 11, 2023
 *      Author: pslocum
 */

#ifndef LMCPOZARRECTANGULARWAVEGUIDE_HH_
#define LMCPOZARRECTANGULARWAVEGUIDE_HH_

#include "LMCField.hh"

#include <vector>

namespace locust
{
 /*!
 @class PozarRectangularWaveguide
 @author P. Slocum
 @brief Derived class to define RectangularWaveguide fields as in Pozar.
 @details
 Available configuration options:
  - "waveguide-central-frequency" -- double [1.63e11] Central frequency in radians
  to be applied in waveguide mode normalizations.

 */

    class PozarRectangularWaveguide: public FieldCore
    {
        public:

	        PozarRectangularWaveguide();
	        virtual ~PozarRectangularWaveguide();

            virtual std::vector<double> TE_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc);
            virtual std::vector<double> TE_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc);
            virtual std::vector<double> TM_E(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc);
            virtual std::vector<double> TM_H(double dimX, double dimY, int m, int n, double xKass, double yKass, double fcyc);

    };


}; /* namespace locust */

#endif /* LMCPOZARRECTANGULARWAVEGUIDE_HH_ */

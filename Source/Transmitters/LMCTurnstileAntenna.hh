/*
 * LMCTurnstileAntenna.hh
 *
 *  Created on: Dec 19, 2019
 *      Author: pslocum
 */

#ifndef LMCTURNSTILEANTENNA_HH_
#define LMCTURNSTILEANTENNA_HH_

#include "LMCTransmitterHardware.hh"

namespace locust
{
 /*!
 @class TurnstileAntenna
 @author P. Slocum
 @brief Derived class describing the orientation and angular dependence of a turnstile antenna transmitter.
 The orientation is such that the antenna sits at (0,0,0) with one dipole moment in the (1,0,0) direction
 and a second dipole in the (0,1,0) direction.  None of these parameters are presently configurable.
 @details
 Available configuration options:
 No input parameters
 */
    class TurnstileAntenna: public TransmitterHardware
    {

        public:
            TurnstileAntenna();
            virtual ~TurnstileAntenna();

            virtual bool Configure( const scarab::param_node& aNode );
            bool Initialize();

            virtual void TxHardwareSayHello();

            virtual double GetPatternFactor(LMCThreeVector pointOfInterest, int antennaNumber);

        private:

            std::vector<LMCThreeVector> fMomentVector;
            std::vector<double> testing;


    };


} /* namespace locust */

#endif /* LMCTURNSTILEANTENNA_HH_ */

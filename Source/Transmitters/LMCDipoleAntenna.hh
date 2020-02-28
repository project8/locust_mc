/*
 * LMCDipoleAntenna.hh
 *
 *  Created on: Dec 19, 2019
 *      Author: pslocum
 */

#ifndef LMCDIPOLEANTENNA_HH_
#define LMCDIPOLEANTENNA_HH_

#include "LMCTransmitterHardware.hh"


namespace locust
{
 /*!
 @class DipoleAntenna
 @author P. Slocum
 @brief Derived class describing the orientation and angular dependence of a dipole antenna transmitter.
 It is assumed that the dipole moment fMomentVector has been oriented such that the electric fields are
 copolar with the receiving antenna elements.  For example, a magnetic dipole moment should typically be
 aligned with the longitudinal axis of an array, and an electric dipole moment is aligned perpendicular
 to an antenna strip and to the longitudinal axis.  If other orientations are selected, then the
 GetPatternFactor() will need to be recalculated.
 @details
 Available configuration options:
 No input parameters
 */
    class DipoleAntenna: public TransmitterHardware
    {

        public:
            DipoleAntenna();
            virtual ~DipoleAntenna();

            virtual bool Configure( const scarab::param_node& aNode );

            virtual void TxHardwareSayHello();

            virtual double GetPatternFactor(LMCThreeVector pointOfInterest, int antennaNumber);

        private:

            LMCThreeVector fMomentVector;
            bool fMagneticDipole;  // either magnetic or electric


    };


} /* namespace locust */

#endif /* LMCDIPOLEANTENNA_HH_ */

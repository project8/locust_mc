/*
 * LMCTransmitterHardware.hh
 *
 *  Created on: Feb. 18, 2020
 *      Author: pslocum
 */

#ifndef LMCTRANSMITTERHARDWARE_HH_
#define LMCTRANSMITTERHARDWARE_HH_

#include "param.hh"
#include "LMCThreeVector.hh"
#include "LMCConst.hh"
#include <vector>

namespace locust
{
 /*!
 @class TransmitterHardware
 @author P. Slocum
 @brief Base class to characterize TransmitterHardware (antenna, orientation in lab) selection
 @details
 Available configuration options:
 - "antenna-x-position": double -- Location of the antenna in the x direction
 - "antenna-y-position": double -- Location of the antenna in the y direction
 - "antenna-z-position": double -- Location of the antenna in the z direction
 No input parameters
 */


    class TransmitterHardware
    {

        public:
            TransmitterHardware();
            virtual ~TransmitterHardware();
            virtual bool Configure( const scarab::param_node& aNode );

            virtual void TxHardwareSayHello();
            
	    int GetNAntennas();

            void SetNAntennas(int aNumber);

            /// Get the positions of the antenna w.r.t the center of the detector
            LMCThreeVector GetAntennaPosition() const;
            
	    /// Set the positions of the antenna w.r.t the center of the detector
            void SetAntennaPosition(const LMCThreeVector &);
            
	    /// Distance from the transmitter to the point provided as argument
	    double GetPropagationDistance(LMCThreeVector pointOfInterest);
            
	    // Factor that describes the radiaiton pattern and depoends on the type of antenna 
            virtual double GetPatternFactor(LMCThreeVector pointOfInterest, int antennaNumber) {return 0.;};

	    // vector pointing from antenna to requested point of interest
            LMCThreeVector ExtractIncidentKVector(LMCThreeVector pointOfInterest);

            virtual double GetDrivePhaseDifference() {return 0.;};

        protected:

	    // Number of antenna, currently set to 1, but could be implemented if needed
            int fNAntennas;
            double fDrivePhaseDifference;

            /// x,y,z positions, could possible be removed if fAntennaPosition is implemented properly
	    double fAntennaPositionX;
            double fAntennaPositionY;
            double fAntennaPositionZ;

	    /// Vector descriving the position of the antenna
            LMCThreeVector fAntennaPosition; // Position of the antenna w.r.t to the center of the array

    };


} /* namespace locust */

#endif

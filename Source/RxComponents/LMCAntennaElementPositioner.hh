/*
 * LMCAntennaElementPositioner.hh
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#ifndef LMCANTENNAELEMENTPOSITIONER_HH_
#define LMCANTENNAELEMENTPOSITIONER_HH_
#include "param.hh"
#include "LMCException.hh"
#include "LMCPatchAntenna.hh"
#include "LMCSlotAntenna.hh"



namespace locust
{
 /*!
 @class LMCAntennaElementPositioner
 @author P. Slocum
 @brief Base class to characterize power combiners
 @details
 Available configuration options:
 No input parameters
 */
    class AntennaElementPositioner
    {

        public:
            AntennaElementPositioner();
            virtual ~AntennaElementPositioner();
            virtual bool Configure( const scarab::param_node& aNode );
            virtual double GetPositionZ(double zShiftArray, int channelIndex, int nChannels,
            		int nSubarrays, int nReceivers, double elementSpacingZ, int receiverIndex);
            double GetTheta(int channelIndex, double dThetaArray);
			virtual void PlaceElement(Receiver &modelElement, double elementRadius, double theta, double zPosition);
            void PlaceAtPosition(Receiver &modelElement, double elementRadius, double theta, double zPosition);
            void PlaceAtConfigured(Receiver &modelElement);
            
        private:
            double fCenterPositionX=0., fCenterPositionY=0., fCenterPositionZ=0.;
            double fPolarizationDirectionX=0., fPolarizationDirectionY=0., fPolarizationDirectionZ=0.;
            double fNormalDirectionX=0., fNormalDirectionY=0., fNormalDirectionZ=0.;
            bool fManualPlacement = false;

};


} /* namespace locust */

#endif

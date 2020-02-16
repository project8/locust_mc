/*
 * LMCTransmitter.hh
 *
 *  Created on: Jan. 24, 2020
 *      Author: pslocum
 */

#ifndef LMCTRANSMITTER_HH_
#define LMCTRANSMITTER_HH_

#include "LMCThreeVector.hh"
#include "LMCChannel.hh"
#include "LMCPatchAntenna.hh"
#include "LMCSlotAntenna.hh"


namespace locust
{
 /*!
 @class Transmitter
 @author P. Slocum
 @brief Base class to characterize Transmitter selection
 @details
 Available configuration options:
 No input parameters
 */


    class Transmitter
    {

        public:
            Transmitter();
            virtual ~Transmitter();
            virtual void TxSayHello();


            virtual double* GetEFieldCoPol(LMCThreeVector pointOfInterest, int channelIndex, int zIndex, double elementSpacing, int nElementsPerStrip, double dt) {};
            virtual LMCThreeVector GetIncidentKVector() {};

            virtual double* SolveKassFields(Receiver* currentElement, double ElementPhi, double tReceiverTime, unsigned tTotalElementIndex) {};
//            virtual double GetEFieldCoPol(Receiver* currentElement, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double ElementPhi) {};
//            virtual double GetEFieldCrossPol(Receiver* currentElement, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double ElementPhi) {};
            virtual void InitializeFieldPoints(std::vector< Channel<Receiver*> > allRxChannels) {};


            virtual bool IsKassiopeia() {return false;};

        private:



};


} /* namespace locust */

#endif /* LMCCHANNEL_HH_ */

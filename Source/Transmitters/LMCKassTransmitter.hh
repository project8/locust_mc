/*
 * LMCKassTransmitter.hh
 *
 *  Created on: Jan 27, 2020
 *      Author: pslocum
 */

#ifndef LMCKASSTRANSMITTER_HH_
#define LMCKASSTRANSMITTER_HH_

#include "LMCTransmitter.hh"
#include "LMCConst.hh"
#include "param.hh"
#include "LMCThreeVector.hh"
#include "LMCLienardWiechert.hh"
#include "LMCGlobalsDeclaration.hh"


namespace locust
{

    /*!
     @class KassTransmitter
     @author P. L. Slocum Jan. 27 2020

     @brief Class to generate tranmistter from plane wave

     @details
     Operates in time space

     Configuration name: "plane-wave"

     Available configuration options:

     */
    class KassTransmitter : public Transmitter
    {
    public:

        KassTransmitter();
        virtual ~KassTransmitter();

        bool Configure( const scarab::param_node& aNode );

        virtual bool IsKassiopeia();

    	double* SolveKassFields(Receiver* currentElement, double ElementPhi, double tReceiverTime, unsigned tTotalElementIndex);
        double GetEFieldCoPol(Receiver* currentElement, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double ElementPhi);
        double GetEFieldCrossPol(Receiver* currentElement, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double ElementPhi);
        double GetAOIFactor(LMCThreeVector IncidentKVector, double ElementPhi);
        void InitializeFieldPoints(std::vector< Channel<Receiver*> > allRxChannels);



    private:

    	LienardWiechert fFieldSolver;


    };


} /* namespace locust */

#endif /* LMCKASSTRANSMITTER_HH_ */

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

namespace locust
{

    /*!
     @class KassTransmitter
     @author P. L. Slocum Jan. 27 2020

     @brief Class to transmit Kassiopeia field solutions

     @details
     Operates in time space

     Configuration name: "transmitter": "kassiopeia"

     Available configuration options:

     */
    class KassTransmitter : public Transmitter
    {
    public:

        KassTransmitter();
        virtual ~KassTransmitter();

        bool Configure( const scarab::param_node& aNode );

        virtual bool IsKassiopeia();

    	std::vector<double> SolveKassFields(LMCThreeVector pointOfInterest, LMCThreeVector coPolDirection, double tReceiverTime, unsigned tTotalElementIndex);
    	void InitializeFieldPoint(LMCThreeVector fieldPoint);

    private:

    	LienardWiechert fFieldSolver;


    };


} /* namespace locust */

#endif /* LMCKASSTRANSMITTER_HH_ */

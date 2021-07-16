/*
 * LMCKassCurrentTransmitter.hh
 *
 *  Created on: June 10, 2021
 *      Author: pslocum
 */

#ifndef LMCKASSCURRENTTRANSMITTER_HH_
#define LMCKASSCURRENTTRANSMITTER_HH_

#include "LMCTransmitter.hh"
#include "LMCConst.hh"
#include "param.hh"
#include "LMCThreeVector.hh"
#include "LMCLienardWiechert.hh"

namespace locust
{

    /*!
     @class KassCurrentTransmitter
     @author P. L. Slocum June 10 2021

     @brief Class to transmit Kassiopeia electron current information.

     @details
     Operates in time space

     Configuration name: "transmitter": "kass-current"

     Available configuration options:

     */
    class KassCurrentTransmitter : public Transmitter
    {
    public:

        KassCurrentTransmitter();
        virtual ~KassCurrentTransmitter();

        bool Configure( const scarab::param_node& aNode );

        virtual bool IsKassiopeia();

        std::vector<double> ExtractParticleXP();


    private:

        double calcOrbitPhase(std::vector<double> tKassParticleXP);
        double quadrantCorrection(double phase, std::vector<double> tKassParticleXP);

    	LienardWiechert fFieldSolver;

    };


} /* namespace locust */

#endif /* LMCKASSCURRENTTRANSMITTER_HH_ */

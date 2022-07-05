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
#include "LMCKassLocustInterface.hh"

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

        std::vector<double> ExtractParticleXP(double TOld, double dt, bool Interpolate, bool Rotate);


    private:

        double calcOrbitPhase(double vx, double vy);
        double calcTheta(double x, double y);
        double quadrantOrbitCorrection(double phase, double vx);
        double quadrantPositionCorrection(double phase, double x);
        int FindNode(double tNew) const;
        double GetGuidingCenterVy();
        double fOrbitPhase;


        kl_interface_ptr_t fInterface;

    };


} /* namespace locust */

#endif /* LMCKASSCURRENTTRANSMITTER_HH_ */

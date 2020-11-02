/*
 * LMCDesignAntenna.hh
 *
 *  Created on: Oct 26, 2020
 *      Author: Arina Telles
 */

#ifndef LMCDESIGNANTENNA_HH_
#define LMCDESIGNANTENNA_HH_

#include "LMCReceiver.hh"

namespace locust
{
 /*!
 @class DesignAntenna
 @author A. Telles
 @brief Derived class describing the an antenna that can have any analytic or externally defined pattern and gain.
 @details
 Intended to be used for antenna design.
 Available configuration options:
 No input parameters
 */
    class DesignAntenna: public Receiver
    {

        public:
            DesignAntenna();
            virtual ~DesignAntenna();

            virtual void RxSayHello();
            virtual double GetPatternFactor(LMCThreeVector incidentKVector, Receiver currentElement);
            virtual void SetFieldPatternExponents(double thetaPatternParameter, double phiPatternParameter);

       	private:
       		double fThetaPatternExponent;
       		double fPhiPatternExponent;

    };


} /* namespace locust */

#endif /* LMCDESIGNANTENNA_HH_ */

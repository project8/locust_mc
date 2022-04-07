/*
 * LMCReceiver.hh
 *
 *  Created on: Mar 1, 2018
 *      Author: nbuzinsky
 */

#ifndef LMCRECEIVER_HH_
#define LMCRECEIVER_HH_

#include "LMCThreeVector.hh"
#include "LMCConst.hh"

namespace locust
{
 /*!
 @class Receiver
 @author N. Buzinsky
 @brief Base class to characterize receiver elements (patches waveguides, etc.)
 @details
 Available configuration options:
 No input parameters
 */
    class Receiver
    {

        public:
            Receiver();
            Receiver(double testvar);
            virtual ~Receiver();

            virtual double GetVoltage() {return 0.;};
            virtual double GetAnalogTimeDelay() {return 0.;};
            virtual LMCThreeVector GetPolarizationDirection();
            virtual void SetPolarizationDirection(const LMCThreeVector &copolDirection);
            virtual LMCThreeVector GetCrossPolarizationDirection();
            virtual void SetCrossPolarizationDirection(const LMCThreeVector &copolDirection);
            virtual LMCThreeVector GetNormalDirection();
            virtual void SetNormalDirection(const LMCThreeVector &normDirection);
            virtual LMCThreeVector GetPosition();
            virtual void SetCenterPosition(const LMCThreeVector &newPosition);
            virtual void RxSayHello();
            virtual double GetPatternFactor(LMCThreeVector incidentKVector, Receiver currentElement);


        private:
            LMCThreeVector copolarizationDirection;
            LMCThreeVector crosspolarizationDirection;
            LMCThreeVector normalDirection;
            LMCThreeVector centerPosition;


};


} /* namespace locust */

#endif /* LMCCHANNEL_HH_ */

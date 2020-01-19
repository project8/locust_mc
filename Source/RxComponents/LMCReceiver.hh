/*
 * LMCReceiver.hh
 *
 *  Created on: Mar 1, 2018
 *      Author: nbuzinsky
 */

#ifndef LMCRECEIVER_HH_
#define LMCRECEIVER_HH_

#include "LMCThreeVector.hh"

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

//            virtual double GetVoltage() = 0;
//            virtual double GetAnalogTimeDelay() = 0;
            LMCThreeVector GetPolarizationDirection();
            void SetPolarizationDirection(const LMCThreeVector &copolDirection);
            LMCThreeVector GetNormalDirection();
            void SetNormalDirection(const LMCThreeVector &normDirection);
            LMCThreeVector GetPosition();
            void SetCenterPosition(const LMCThreeVector &newPosition);


        private:
            LMCThreeVector copolarizationDirection;
            LMCThreeVector normalDirection;
            LMCThreeVector centerPosition;


};


} /* namespace locust */

#endif /* LMCCHANNEL_HH_ */

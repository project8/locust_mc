/*
 * LMCPatchAntenna.hh
 *
 *  Created on: Mar 1, 2018
 *      Author: nbuzinsky
 */

#ifndef LMCPATCHANTENNA_HH_
#define LMCPATCHANTENNA_HH_

#include <boost/math/interpolators/cubic_b_spline.hpp>
#include "LMCThreeVector.hh"
#include "LMCConst.hh"
#include "LMCReceiver.hh"

namespace locust
{
 /*!
 @class PatchAntenna
 @author N. Buzinsky
 @brief Derived class describing the behaviours of the patch antenna voltage, as a function of the E field
 @details
 Available configuration options:
 No input parameters
 */
    class PatchAntenna: public Receiver
    {

        public:
            PatchAntenna();
            PatchAntenna(const LMCThreeVector &patchPosition);
            PatchAntenna(const LMCThreeVector &patchPosition, const double &timeDelay);
            virtual ~PatchAntenna();

            virtual double GetVoltage();
            virtual double GetAnalogTimeDelay();


            void SetIncidentElectricField(const LMCThreeVector &incomingElectricField);
            void SetIncidentMagneticField(const LMCThreeVector &incomingMagneticField);
            void SetInstantaneousFrequency(const double &dopplerFrequency);

            LMCThreeVector GetPosition();

            int GetPreviousRetardedIndex();
            double  GetPreviousRetardedTime();

            void SetPreviousRetardedIndex(const int& index);
            void SetPreviousRetardedTime(const double &time);

            void SetCenterPosition(const LMCThreeVector &newPosition);
            void SetPolarizationDirection(const LMCThreeVector &copolDirection);



        private:
            double GetAntennaFactor();
            double GetGainFactor(); 
            double GetCopolarizationFactor();

            LMCThreeVector copolarizationDirection;
            LMCThreeVector normalDirection;

            LMCThreeVector centerPosition;

            LMCThreeVector incidentElectricField;
            LMCThreeVector incidentMagneticField;

            int previousRetardedIndex;
            double previousRetardedTime;

            double instantaneousFrequency;

            double timeDelay;

            const std::vector<double> antennaFactor = {877.12321,626.12312,608.12312,620.77,905.37};
            const double lowerBoundFrequency = 25.1e9;
            const double frequencySpacingSpline = 0.1e9;

            const std::vector<double> gain= {4.010465406, 3.858038578, 3.399096589, 2.73347503, 2.009885519, 1.364161623, 0.869458575, 0.533070658, 0.324008154, 0.201512165};
            const double lowerBoundAngle = 0;
            const double angularSpacingSpline = 10 * LMCConst::Pi() / 180.; //Angular spacing in radians

            boost::math::cubic_b_spline<double> antennaFactorSpline;
            boost::math::cubic_b_spline<double> gainSpline;

    };


} /* namespace locust */

#endif /* LMCPATCHANTENNA_HH_ */

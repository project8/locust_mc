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
            void SetNormalDirection(const LMCThreeVector &normDirection);


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

            const std::vector<double> antennaFactor = {577.051, 571.497, 565.992, 560.537, 555.136, 549.789, 544.499, 539.269, 534.100, 528.995, 523.957, 518.988, 514.091, 509.268, 504.523, 499.857, 495.275, 490.779, 486.372, 482.058, 477.839, 473.719, 469.702, 465.790, 461.987, 458.296, 454.722, 451.268, 447.936, 444.732, 441.658, 438.718, 435.916, 433.255, 430.739, 428.371, 426.155, 424.093, 422.189, 420.446, 418.868, 417.456, 416.214, 415.144, 414.248, 413.529, 412.988, 412.627, 412.448, 412.451, 412.638, 413.010, 413.568, 414.310, 415.238, 416.352, 417.650, 419.133, 420.798, 422.646, 424.675, 426.884, 429.270, 431.831, 434.567, 437.474, 440.550, 443.792, 447.199, 450.768, 454.495, 458.378, 462.413, 466.599, 470.932, 475.409, 480.027, 484.784, 489.675, 494.698, 499.850, 505.129, 510.530, 516.052, 521.691, 527.445, 533.311, 539.286, 545.368, 551.553, 557.841, 564.227, 570.710, 577.287, 583.956, 590.715, 597.561, 604.494, 611.509, 618.607};

            const double lowerBoundFrequency = 25.1e9;
            const double frequencySpacingSpline = 1.e9 /99.;

            const std::vector<double> gain= { 413.2552, 413.4029, 413.7279, 414.2307, 414.9122, 415.7736, 416.8165, 418.0426, 419.4540, 421.0532, 422.8427, 424.8257, 427.0053, 429.3850, 431.9687, 434.7602, 437.7639, 440.9843, 444.4262, 448.0945, 451.9944, 456.1315, 460.5115, 465.1403, 470.0241, 475.1694, 480.5829, 486.2715, 492.2426, 498.5036, 505.0624, 511.9270, 519.1060, 526.6080, 534.4421, 542.6178, 551.1447, 560.0330, 569.2932, 578.9360, 588.9725, 599.4145, 610.2737, 621.5624, 633.2932, 645.4791, 658.1333, 671.2694, 684.9012, 699.0430, 713.7091, 728.9140, 744.6726, 760.9996, 777.9100, 795.4189, 813.5410, 832.2914, 851.6845, 871.7349, 892.4566, 913.8632, 935.9679, 958.7832, 982.3209, 1006.591, 1031.605, 1057.372, 1083.897, 1111.189, 1139.251, 1168.086, 1197.697, 1228.081, 1259.237, 1291.159, 1323.840, 1357.270, 1391.437, 1426.325, 1461.917, 1498.192, 1535.126, 1572.693, 1610.862, 1649.602, 1688.878, 1728.650, 1768.877, 1809.517, 1850.522};

            const double lowerBoundAngle = 0;
            const double angularSpacingSpline = 1. * LMCConst::Pi() / 180.; //Angular spacing in radians

            boost::math::cubic_b_spline<double> antennaFactorSpline;
            boost::math::cubic_b_spline<double> gainSpline;

    };


} /* namespace locust */

#endif /* LMCPATCHANTENNA_HH_ */

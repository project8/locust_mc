
#ifndef LMCPOWERCOMBINER_HH_
#define LMCPOWERCOMBINER_HH_

#include "LMCConst.hh"
#include <vector>

namespace locust
{
 /*!
 @class PowerCombiner
 @author P. Slocum
 @brief Class to describe the power combining in the patch array.
 @details
 Available configuration options:
 No input parameters
 */
    class PowerCombiner
    {

        public:
            PowerCombiner();

            virtual ~PowerCombiner();

            double GetCorporateVoltageDamping();
            double GetSeriesPhaseDelay(unsigned z_index, double DopplerFrequency, double PatchSpacing);
            double GetSeriesVoltageDamping(unsigned z_index);
            double GetCenterFedPhaseDelay(int NPatchesPerStrip, unsigned z_index, double DopplerFrequency, double PatchSpacing);
            double GetOneQuarterVoltageDamping(int NPatchesPerStrip, unsigned z_index);
            double GetSevenEighthsVoltageDamping(int NPatchesPerStrip, unsigned z_index);
            double GetNineSixteenthsVoltageDamping(int NPatchesPerStrip, unsigned z_index);
            double GetVoltageDividerWeight(double RJunction, double R0, double Rground, int NPatchesPerStrip, unsigned z_index);
            std::vector<double> GetResistances(double RJunction, double R0, double RGround, int NPAIRS);
            std::vector<double> GetPartialGains(double RJunction, double R0, double RGround, int NPAIRS);
            double GetParallelResistance(std::vector<double> R, int NRESISTORS, int resistorindex);


        private:

    };


} /* namespace locust */

#endif /* LMCPOWERCOMBINER_HH_ */


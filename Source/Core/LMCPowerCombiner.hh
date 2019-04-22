
#ifndef LMCPOWERCOMBINER_HH_
#define LMCPOWERCOMBINER_HH_

#include "LMCConst.hh"

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

        private:

    };


} /* namespace locust */

#endif /* LMCPOWERCOMBINER_HH_ */


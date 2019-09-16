
#ifndef LMCPOWERCOMBINER_HH_
#define LMCPOWERCOMBINER_HH_

#include "LMCSignal.hh"
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

            bool AddOneVoltageToStripSum(Signal* aSignal, double VoltageFIRSample, double phi_LO, unsigned z_index, unsigned sampleIndex);
            double GetCorporateVoltageDamping();
            double GetSeriesPhaseDelay(unsigned z_index, double DopplerFrequency, double PatchSpacing);
            double SetSeriesVoltageDamping();
            double GetCenterFedPhaseDelay(int NPatchesPerStrip, unsigned z_index, double DopplerFrequency, double PatchSpacing);
            double GetOneQuarterVoltageDamping(int NPatchesPerStrip, unsigned z_index);
            double GetSevenEighthsVoltageDamping(int NPatchesPerStrip, unsigned z_index);
            double GetNineSixteenthsVoltageDamping(int NPatchesPerStrip, unsigned z_index);
            double GetSeriesVoltageDamping(unsigned z_index);
            double GetVoltageDividerWeight(double RJunction, double R0, double Rground, int NPatchesPerStrip, unsigned z_index);
            void SetVoltageDampingFactors(int aPowerCombiner, int aPatchesPerStrip);
            void SetSMatrixParameters(int powerCombiner, int aPatchesPerStrip);
            void SetNPatchesPerStrip(int aPatchesPerStrip);
            void SetJunctionLoss(double aJunctionLoss);
            void SetPatchLoss(double aPatchLoss);
            void SetAmplifierLoss(double aAmplifierLoss);
            void SetEndPatchLoss(double aEndPatchLoss);



        private:
            void SetCenterFedDampingFactors();
            void SetSeriesFedDampingFactors();
            void SetVoltageDividerDampingFactors();
            std::vector<double> GetResistances(double RJunction, double R0, double RGround, int NPAIRS);
            std::vector<double> GetPartialGains(double RJunction, double R0, double RGround, int NPAIRS);
            double GetParallelResistance(std::vector<double> R, int NRESISTORS, int resistorindex);
            int nPatchesPerStrip;
            double junctionLoss;
            double patchLoss;
            double amplifierLoss;
            double endPatchLoss;
            double junctionResistance;
      	    std::vector<double> dampingFactors;


    };


} /* namespace locust */

#endif /* LMCPOWERCOMBINER_HH_ */


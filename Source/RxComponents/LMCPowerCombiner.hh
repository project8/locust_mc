
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
 S-matrix and voltage damping factors for the array are defined with functions
 SetSMatrixParameters(int powerCombiner, int aPatchesPerStrip) and
 SetVoltageDampingFactors(int aPowerCombiner, int aPatchesPerStrip) .

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
            void SetVoltageDampingFactors(int aPowerCombiner, int aPatchesPerStrip);
            void SetSMatrixParameters(int powerCombiner, int aPatchesPerStrip);
            void SetNPatchesPerStrip(int aPatchesPerStrip);
            void SetJunctionLoss(double aJunctionLoss);
            void SetPatchLoss(double aPatchLoss);
            void SetAmplifierLoss(double aAmplifierLoss);
            void SetEndPatchLoss(double aEndPatchLoss);



        private:
            double GetSeriesPhaseDelay(unsigned z_index, double DopplerFrequency, double PatchSpacing);
            double GetCenterFedPhaseDelay(int NPatchesPerStrip, unsigned z_index, double DopplerFrequency, double PatchSpacing);
            void SetCenterFedDampingFactors();
            void SetSeriesFedDampingFactors();
            void SetVoltageDividerDampingFactors();
            std::vector<double> GetResistances(double RJunction, double R0, double RGround, int NPAIRS);
            std::vector<double> GetPartialGains(double RJunction, double R0, double RGround, int NPAIRS);
            double GetVoltageDividerWeight(double RJunction, double R0, double Rground, int NPatchesPerStrip, unsigned z_index);
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


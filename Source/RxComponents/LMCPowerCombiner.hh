
#ifndef LMCPOWERCOMBINER_HH_
#define LMCPOWERCOMBINER_HH_

#include "LMCException.hh"
#include "param.hh"
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
 SetSMatrixParameters(int aPatchesPerStrip) and
 SetVoltageDampingFactors(int aPatchesPerStrip) .

 @details
 Available configuration options:
 "power-combining-feed" is an integer to select the appropriate power combining configuration.
 No input parameters
 */
    class PowerCombiner
    {

        public:
            PowerCombiner();

            virtual ~PowerCombiner();
            bool Configure( const scarab::param_node& aNode);

            bool AddOneVoltageToStripSum(Signal* aSignal, double VoltageFIRSample, double phi_LO, unsigned z_index, unsigned sampleIndex);
            bool SetVoltageDampingFactors(int aPatchesPerStrip);
            bool SetSMatrixParameters(int aPatchesPerStrip);
            void SetNPatchesPerStrip(int aPatchesPerStrip);
            void SetJunctionLoss(double aJunctionLoss);
            void SetPatchLoss(double aPatchLoss);
            void SetAmplifierLoss(double aAmplifierLoss);
            void SetEndPatchLoss(double aEndPatchLoss);
            bool SetPowerCombiner( std::string feed );
            int GetPowerCombiner();



        private:
            double GetSeriesPhaseDelay(unsigned z_index, double DopplerFrequency, double PatchSpacing);
            double GetCenterFedPhaseDelay(unsigned z_index, double DopplerFrequency, double PatchSpacing);
            bool SetCenterFedDampingFactors();
            bool SetSeriesFedDampingFactors();
            bool SetVoltageDividerDampingFactors();
            bool SetSmatrix10patchDampingFactors();
            std::vector<double> GetResistances(double RJunction, double R0, double RGround, int NPAIRS);
            std::vector<double> GetPartialGains(double RJunction, double R0, double RGround, int NPAIRS);
            double GetVoltageDividerWeight(double RJunction, double R0, double Rground, unsigned z_index);
            double GetParallelResistance(std::vector<double> R, int NRESISTORS, int resistorindex);
            int fpowerCombiner;
            int fnPatchesPerStrip;
            double fjunctionLoss;
            double fpatchLoss;
            double famplifierLoss;
            double fendPatchLoss;
            double fjunctionResistance;
      	    std::vector<double> fdampingFactors;

      	    std::vector<double> fsMatrix10patch = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

    };


} /* namespace locust */

#endif /* LMCPOWERCOMBINER_HH_ */


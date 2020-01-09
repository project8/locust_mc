
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
     	if (feed == "corporate") fpowerCombiner = 0;  // default
    	else if (feed == "series") fpowerCombiner = 1;
        else if (feed == "one-quarter") fpowerCombiner = 2;
        else if (feed == "seven-eighths") fpowerCombiner = 3;
        else if (feed == "nine-sixteenths") fpowerCombiner = 4;
        else if (feed == "voltage-divider") fpowerCombiner = 5;
        else if (feed == "s-matrix") fpowerCombiner = 6;
        else if (feed == "single-patch") fpowerCombiner = 7;
        else fpowerCombiner = 0;  // default

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
            bool SetSmatrixDampingFactors();
            std::vector<double> GetResistances(double RJunction, double R0, double RGround, int NPAIRS);
            std::vector<double> GetPartialGains(double RJunction, double R0, double RGround, int NPAIRS);
            double GetVoltageDividerWeight(double RJunction, double R0, double Rground, unsigned z_index);
            double GetParallelResistance(std::vector<double> R, int NRESISTORS, int resistorindex);
            std::vector<double> GetSmatrixElements();
            bool SetTransmissionCoefficients();
            int fpowerCombiner;
            int fnPatchesPerStrip;
            double fjunctionLoss;
            double fpatchLoss;
            double famplifierLoss;
            double fendPatchLoss;
            double fjunctionResistance;
      	    std::vector<double> fdampingFactors;

      	    // Uniform taper S-matrices from HFSS:
      /*
      	    std::vector<double> fsMatrix2patch = {0.2, 0.64, 0.64};
            std::vector<double> fsMatrix4patch = {0.09, 0.47, 0.47, 0.47, 0.47};
            std::vector<double> fsMatrix6patch = {0.03, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38};
      */
            // end Uniform taper S-matrices.

            // 7/8 power combiner S-matrices (traveling wave configuration) from HFSS:
            std::vector<double> fsMatrix2patch = {0.12, 0.43, 0.43};
            std::vector<double> fsMatrix4patch = {0.12, 0.3, 0.43, 0.43, 0.3};
            std::vector<double> fsMatrix6patch = {0.12, 0.24, 0.3, 0.43, 0.43, 0.3, 0.24};
            std::vector<double> fsMatrix8patch = {0.12, 0.17, 0.24, 0.3, 0.43, 0.43, 0.3, 0.24, 0.17};
            // end 7/8 combiner S-matrices.

            std::vector<double> ftransmissionCoefficients;
    };


} /* namespace locust */

#endif /* LMCPOWERCOMBINER_HH_ */


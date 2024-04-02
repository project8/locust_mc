/*
 * LMCPowerCombiner.hh
 *
 *  Created on: Feb 25, 2020
 *      Author: pslocum
 */

#ifndef LMCPOWERCOMBINER_HH_
#define LMCPOWERCOMBINER_HH_
#include "param.hh"
#include "LMCException.hh"
#include "LMCConst.hh"
#include "LMCSignal.hh"
#include "LMCPatchAntenna.hh"
#include "LMCSlotAntenna.hh"
#include <vector>


namespace locust
{
 /*!
 @class LMCPowerCombiner
 @author P. Slocum
 @brief Base class to characterize power combiners
 @details
 Available configuration options:
      - "waveguide-short":  bool [true] optional presence/absence of reflecting short in waveguide.
      - "voltage-check": bool [false] optional print signal voltages to terminal during simulation.
      This can be useful for e.g. checking suitability of digitizer range or other debugging.
 */
    class PowerCombiner
    {

        public:
            PowerCombiner();
            virtual ~PowerCombiner();
            int GetNElementsPerStrip();
            void SetNElementsPerStrip( int aNumberOfElements );
            virtual bool Configure( const scarab::param_node& aNode );
        	virtual bool SetVoltageDampingFactors() {return true;};
        	virtual bool SetSMatrixParameters() {return true;};
        	virtual bool IsSinglePatch();
            virtual Receiver* ChooseElement();
        	bool AddOneVoltageToStripSum(Signal* aSignal, double excitationAmplitude, double phi_LO, unsigned z_index, unsigned sampleIndex);
        	virtual bool AddOneModeToCavityProbe(int l, int m, int n, Signal* aSignal, std::vector<double> particleXP, double excitationAmplitude, double EFieldAtProbe, std::vector<double> dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, int channelIndex, bool initParticle) {return true;};
		virtual bool AddOneSampleToRollingAvg(int bTE, int l, int m, int n, double excitationAmplitude, unsigned sampleIndex) {return true;};
        	virtual void SayHello();
        	virtual void Initialize() {};



            double GetJunctionLoss();
            void SetJunctionLoss( double aJunctionLoss );
            double GetPatchLoss();
            void SetPatchLoss( double aPatchLoss );
            double GetAmplifierLoss();
            void SetAmplifierLoss( double aAmplifierLoss );
            double GetEndPatchLoss();
            void SetEndPatchLoss( double aEndPatchLoss );
            double GetJunctionResistance();
            void SetJunctionResistance( double aJunctionResistance );
            double GetDampingFactor( int z_index );
            void SetDampingFactor (int z_index, double aDampingFactor );
            bool GetVoltageCheck();
            void SetNCavityModes( int aNumberOfModes );
            int GetNCavityModes();
            void SetNChannels( int aNumberOfChannels );
            int GetNChannels();
            void SetOutputPath( std::string aPath );
            std::string GetOutputPath();
            double GetVoltagePhase();
            void SetVoltagePhase( double aPhase );
            virtual bool SizeNChannels(int aNumberOfChannels) {return true;};

            bool GetWaveguideShortIsPresent();
            void SetWaveguideShortIsPresent ( bool aValue );



        private:
            int fnElementsPerStrip;
      	    std::vector<double> fdampingFactors;
            double fjunctionLoss;
            double fpatchLoss;
            double famplifierLoss;
            double fendPatchLoss;
            double fjunctionResistance;
            bool fvoltageCheck;
            int fNCavityModes;
            int fNChannels;
            std::string fOutputPath;
            bool fWaveguideShortIsPresent;



};


} /* namespace locust */

#endif

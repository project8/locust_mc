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
 No input parameters
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
        	virtual bool AddOneModeToCavityProbe(Signal* aSignal, std::vector<double> particleXP, double excitationAmplitude, double EFieldAtProbe, std::vector<double> dopplerFrequency, double dt, double phi_LO, double totalScalingFactor, unsigned sampleIndex, bool initParticle) {return true;};
        	virtual bool AddOneSampleToRollingAvg(int l, int m, int n, double excitationAmplitude, unsigned sampleIndex) {return true;};
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
            double GetCavityProbeZ();
            void SetCavityProbeZ ( double aZ );
            double GetCavityProbeRFrac();
            void SetCavityProbeRFrac ( double aFraction );
            double GetVoltagePhase();
            void SetVoltagePhase( double aPhase );

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
            double fCavityProbeZ;
            double fCavityProbeRFrac;
            bool fWaveguideShortIsPresent;



};


} /* namespace locust */

#endif

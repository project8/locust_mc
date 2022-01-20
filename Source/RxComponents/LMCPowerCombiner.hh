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
        	bool AddOneVoltageToStripSum(Signal* aSignal, double VoltageFIRSample, double phi_LO, unsigned z_index, unsigned sampleIndex);
        	bool AddOneModeToCavityProbe(Signal* aSignal, double VoltageFIRSample, double phi_LO, double totalScalingFactor, double cavityProbeImpedance, unsigned sampleIndex);
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
            std::vector<double> GetCavityProbeZ();
            void SetCavityProbeZ ( std::vector<double> aVector );
            std::vector<double> GetCavityProbeTheta();
            void SetCavityProbeTheta ( std::vector<double> aVector );
            int GetNCavityProbes();
            void SetNCavityProbes( int aNumberOfProbes );
            void SetNCavityModes( int aNumberOfModes );
            double GetCavityProbeInductance();
            void SetCavityProbeInductance( double anInductance );
            bool SetCavityProbeLocations(int nCavityProbes, double cavityLength);
        	bool AddOneSampleToRollingAvg(int l, int m, int n, double VoltageFIRSample, double totalScalingFactor, unsigned sampleIndex);


        private:
            int fnElementsPerStrip;
      	    std::vector<double> fdampingFactors;
            double fjunctionLoss;
            double fpatchLoss;
            double famplifierLoss;
            double fendPatchLoss;
            double fjunctionResistance;
            bool fvoltageCheck;
            int fnCavityProbes;
            int fNCavityModes;
            double fCavityProbeInductance;
            std::vector<double> fCavityProbeZ;
            std::vector<double> fCavityProbeTheta;
            std::vector<std::vector<std::vector<double>>> fRollingAvg;
            std::vector<std::vector<std::vector<int>>> fCounter;


};


} /* namespace locust */

#endif

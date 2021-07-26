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
        	virtual bool SetVoltageDampingFactors() {};
        	virtual bool SetSMatrixParameters() {};
        	virtual bool IsSinglePatch();
            virtual Receiver* ChooseElement();
        	bool AddOneVoltageToStripSum(Signal* aSignal, double VoltageFIRSample, double phi_LO, unsigned z_index, unsigned sampleIndex);
        	bool AddOneModeToCavityProbe(Signal* aSignal, double VoltageFIRSample, double phi_LO, double modeAmplitudeFactor, double cavityProbeImpedance, unsigned sampleIndex);
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
            double GetCavityProbeImpedance();
            void SetCavityProbeImpedance( double anImpedance );
            bool SetCavityProbeLocations(int nCavityProbes, double cavityLength);


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
            double fCavityProbeImpedance;
            std::vector<double> fCavityProbeZ;
            std::vector<double> fCavityProbeTheta;

};


} /* namespace locust */

#endif

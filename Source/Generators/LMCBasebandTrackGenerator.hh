/*
 * LMCBasebandTrackGenerator.hh
 *
 *  Created on: Jan 14 2015
 *      Author: plslocum after nsoblath
 */

#ifndef LMCBASEBANDTRACKGENERATOR_HH_
#define LMCBASEBANDTRACKGENERATOR_HH_

#include "../Core/LMCGenerator.hh"
#include "../Core/LMCRunLengthCalculator.hh"
#define PI 3.1415926 // pls:  this should go into ../Core/LMCConstants.hh .


namespace locust
{

    /*!
     @class BasebandTrackGenerator
     @author P. L. Slocum after N. S. Oblath

     @brief Add Sine Wave to the signal.  Choose starting electron 
       energy and then lose energy to Larmor power.

     @details
     Can operate in time or frequency space

     Configuration name: "baseband-track"

     Available configuration options:
     - "electron-energy": double -- Energy of electron (keV).
     - "total-lo-freqs": double -- Sum of LO frequencies (GHz).
     - "domain": string -- Determines whether the sinusoidal test signal is generated in the time 
            or frequency domain
    
     Available options: "time" and "freq" [default]

    */
    class BasebandTrackGenerator : public Generator
    {
        public:
            BasebandTrackGenerator( const std::string& aName = "baseband-track" );
            virtual ~BasebandTrackGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            double GetElectronEnergy() const;
            void SetElectronEnergy( double aElectronEnergy );

            double GetTotalLOFreqs() const;
            void SetTotalLOFreqs( double aTotalLOFreqs );

            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );

        private:
            bool DoGenerate( Signal* aSignal ) const;

            bool DoGenerateTime( Signal* aSignal ) const;
            bool DoGenerateFreq( Signal* aSignal ) const;

            bool (BasebandTrackGenerator::*fDoGenerateFunc)( Signal* aSignal ) const;
            double CalculateLarmorPower( double gamma ) const;
            double CalculateGamma( double KineticEnergy ) const;
            double CalculateCyclotronFrequency( double Gamma ) const;
            double CalculateBasebandFrequency( double CyclotronFrequency ) const;

            double fElectronEnergy;
            double fTotalLOFreqs;

            
    };

} /* namespace locust */

#endif /* LMCBasebandTrackGENERATOR_HH_ */


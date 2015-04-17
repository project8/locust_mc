/*
 * LMCTRAPPEDELECTRONGENERATOR.hh
 *
 *  Created on: Mar 6, 2015
 *      Author: plslocum after nsoblath
 */

#ifndef LMCTRAPPEDELECTRONGENERATOR_HH_
#define LMCTRAPPEDELECTRONGENERATOR_HH_

#include "../Core/LMCGenerator.hh"
#include "../Core/LMCRunLengthCalculator.hh"
#define LO_FREQUENCY 26.385015e9 // Hz
#define PI 3.1415926
#define Z1 -2.5 // z-offset (cm) of lower coil.
#define Z2 2.5  // z-offset (cm) of upper coil.
#define RCOIL 2.0 // radius (cm) of both coils.
#define CENTER_TO_ANTENNA 11.0 // distance from antenna to center of trap (cm).
#define CENTER_TO_SHORT 11.0 // distance from short to center of trap (cm).
#define C 2.99792458e10 // speed of light in cm/s.
#define MU0 1.256637e-4 // magnetic constant in T*cm/A.
#define CURRENT 2.0 // bathtub coil current in amps.


namespace locust
{

    /*!
     @class TrappedElectronGenerator
     @author P. L. Slocum after N. S. Oblath

     @brief Add Gaussian-distributed noise to the signal

     @details
     Operates in time space

     Configuration name: "trapped-electron"

     Available configuration options:

    */
    class TrappedElectronGenerator : public Generator
    {
        public:
            TrappedElectronGenerator( const std::string& aName = "trapped-electron" );
            virtual ~TrappedElectronGenerator();

            bool Configure( const ParamNode* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            Signal::State GetDomain() const;
            void SetDomain( Signal::State aDomain );

        private:
            bool DoGenerate( Signal* aSignal ) const;

            bool DoGenerateTime( Signal* aSignal ) const;
            bool DoGenerateFreq( Signal* aSignal ) const;
            bool (TrappedElectronGenerator::*fDoGenerateFunc)( Signal* aSignal ) const;

            double GetBMag(double x0, double y0, double z0) const;
            double GetBz( double x0, double y0, double z0 ) const;
            double GetBx( double x0, double y0, double z0 ) const;
            double GetBy( double x0, double y0, double z0 ) const;
            double GetSpeed(double KineticEnergy) const;
            double GetKineticEnergy(double Velocity) const;
            double CalculateGamma( double KineticEnergy ) const;
            double *StepElectron(double *OldPosition, double KineticEnergy, double mu, double dt) const;
            double *StartElectron(double KineticEnergy, double PitchAngle) const;
            bool StopElectron(double *position) const;
            double GetMu(double *position, double KineticEnergy, double PitchAngle) const;
            double CalculateLarmorPower( double gamma, double B) const;  
            double CalculateCyclotronFrequency( double Gamma , double *position) const;
            double CalculateBasebandFrequency( double CyclotronFrequency ) const;
            double GetVoltagePhase( double Freq, double dt ) const;
            double GetVoltagePhaseFromShort( double Freq, double dt ) const;
            double GetCyclotronFreqAntennaFrame( double RFFreq, double Vparallel) const;

    };

} /* namespace locust */

#endif /* LMCTRAPPEDELECTRONGENERATOR_HH_ */

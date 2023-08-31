/*
 * LMCDampedHarmonicOscillator.hh
 *
 *  Created on: Jul 6, 2022
 *      Author: pslocum
 *
This class calculates an analytic Green's function from a model of a damped harmonic
oscillator.  To use the model in the context of the LMCCavitySignalGenerator, the
following default syntax can be used:

LocustSim -c config/LocustCavityTemplate.json

To print out the resulting FIR to a text file "output/FIR.txt", use this parameter:

"cavity-signal.print-fir-debug"=true
 *
 *
 *
 */

#ifndef LMCDAMPEDHARMONICOSCILLATOR_HH_
#define LMCDAMPEDHARMONICOSCILLATOR_HH_
#include "LMCAnalyticResponseFunction.hh"
#include "LMCConst.hh"
#include "param.hh"
#include "LMCException.hh"
#include <math.h>
#include "LMCFIRFileHandler.hh"
#include "LMCTFFileHandler.hh"
#include "LMCSignal.hh"


namespace locust
{
    class DampedHarmonicOscillator : public AnalyticResponseFunction
    {

        public:
    		DampedHarmonicOscillator();
    		virtual ~DampedHarmonicOscillator();
    		virtual bool Configure( const scarab::param_node& aNode );
            virtual bool GenerateGreensFunction(int bTE, int l, int m, int n);
            bool Initialize();
            std::pair<double,double> GreensFunction(int bTE, int l, int m, int n, double t);
            double ExpDecayTerm(int bTE, int l, int m, int n, double t);
            virtual void SetCavityQ( int bTE, int l, int m, int n, double aQ );
            virtual double GetCavityQ(int bTE, int l, int m, int n);
            virtual void SetCavityFrequency( int bTE, int l, int m, int n, double aFrequency );
            virtual double GetCavityFrequency(int bTE, int l, int m, int n);
            virtual void SetDHOTimeResolution(int bTE, int l, int m, int n, double aTimeResolution );
            virtual double GetDHOTimeResolution(int bTE, int l, int m, int n);
            virtual void SetDHOThresholdFactor(int bTE, int l, int m, int n, double aThresholdFactor );
            virtual double GetDHOThresholdFactor(int bTE, int l, int m, int n);
            double NormFactor(int bTE, int l, int m, int n, double aDriveFrequency);
            bool PopulateCalibrationSignal(int bTE, int l, int m, int n, Signal* aSignal, int N0, double aDriveFrequency);
            std::deque<double> SignalToDeque(Signal* aSignal);
	    std::deque<double> SignalToDequeArray(int bTE, int l, int m, int n, Signal* aSignal);
	    //virtual int GetNModes();


        private:
            int fNModes;
            std::vector<std::vector<std::vector<std::vector<double>>>> fCavityFrequency; // Hz
            std::vector<std::vector<std::vector<std::vector<double>>>> fCavityOmega; // radians/s
            std::vector<std::vector<std::vector<std::vector<double>>>> fCavityQ;
            std::vector<std::vector<std::vector<std::vector<double>>>> fCavityDampingFactor;
            std::vector<std::vector<std::vector<std::vector<double>>>> fBFactor; // harmonic oscillator parameter
            std::vector<std::vector<std::vector<std::vector<double>>>> fCavityOmegaPrime;  // damped resonant frequency
            int fMaxNBins;
            std::vector<std::vector<std::vector<std::vector<double>>>> fTimeResolution;
            std::vector<std::vector<std::vector<std::vector<double>>>> fThresholdFactor;
            std::vector<std::vector<std::vector<std::vector<double>>>> fHannekePowerFactor;
            TFReceiverHandler* fTFReceiverHandler;


    };


} /* namespace locust */

#endif

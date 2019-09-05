/*
 * LMCButterworthLPFGenerator.cc
 *
 *  Created on: Sept 9, 2016
 *      Author: plslocum
 */

#include "LMCButterworthLPFGenerator.hh"

#include "logger.hh"

using std::string;

namespace locust
{
    LOGGER( lmclog, "ButterworthLPFGenerator" );

    MT_REGISTER_GENERATOR(ButterworthLPFGenerator, "butterworth-lpf");

    ButterworthLPFGenerator::ButterworthLPFGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &ButterworthLPFGenerator::DoGenerateTime ),
		fFIRa( 0. ),
		fFIRb( 0. ),
		fN( 0. ),
		fM( 0. )
    {
        fRequiredSignalState = Signal::kTime;
    }

    ButterworthLPFGenerator::~ButterworthLPFGenerator()
    {
    }

    bool ButterworthLPFGenerator::SetCoefficients()
    {

    	fN = 8;
    	fM = 8;
    	fFIRa.resize(fN+1);
    	fFIRb.resize(fN+1);

        // raw coefficients from numerator of G(z) in Mathematica, starting with Butterworth polynomials.
        fFIRa[0] = 12.675;
        fFIRa[1] = 101.4;
        fFIRa[2] = 354.9;
        fFIRa[3] = 709.8;
        fFIRa[4] = 887.25;
        fFIRa[5] = 709.8;
        fFIRa[6] = 354.9;
        fFIRa[7] = 101.4;
        fFIRa[8] = 12.675;

        // raw coefficients from denominator of G(z) in Mathematica, starting with Butterworth polynomials.
        fFIRb[0] = 160.698;
        fFIRb[1] = -511.67;
        fFIRb[2] = -832.628;
        fFIRb[3] = -837.901;
        fFIRb[4] = -561.456;
        fFIRb[5] = -252.64;
        fFIRb[6] = -73.9998;
        fFIRb[7] = -12.8136;
        fFIRb[8] = -1.;


        double norm = fFIRb[0];
        for (int i=0; i<fN+1; i++)
          {
          fFIRa[i] /= norm;
          fFIRb[i] /= norm;
          }

    	return true;
    }

    bool ButterworthLPFGenerator::Configure( const scarab::param_node& aParam )
    {

        SetCoefficients();

    	return true;
    }

    void ButterworthLPFGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    bool ButterworthLPFGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool ButterworthLPFGenerator::DoGenerateTime( Signal* aSignal )
    {
        const unsigned nchannels = fNChannels;

        for (int ch=0; ch<nchannels; ch++)
        {
        for( unsigned IQindex = 0; IQindex < 2; IQindex++ )
        {

        double* FilteredMagnitude = new double[aSignal->TimeSize()]; // temporary signal.
        double* RawMagnitude = new double[aSignal->TimeSize()];  // incoming signal.

        for( unsigned index = 0; index < aSignal->TimeSize(); ++index ) // initialize
          {
          FilteredMagnitude[index] = 0.;
          RawMagnitude[index] = aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][IQindex];
          }

        for( unsigned index = fN; index < aSignal->TimeSize(); ++index )
          {

          double FIR_component = 0.;
          double IIR_component = 0.;

          for (int i=0; i<fN+1; i++)
            {
            FIR_component += fFIRa[i]*RawMagnitude[index-i];
            }
          for (int i=1; i<fM+1; i++)
            {
            IIR_component += fFIRb[i]*FilteredMagnitude[index-i];
            }

          FilteredMagnitude[index] = FIR_component + IIR_component;

          }

        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )  // apply filter
          {
          aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index][IQindex] = FilteredMagnitude[index];
          }

        delete FilteredMagnitude;
        delete RawMagnitude;

        } // IQindex
        } // nchannels



        return true;
    }

    bool ButterworthLPFGenerator::DoGenerateFreq( Signal* aSignal )
    {

        return true;
    }

} /* namespace locust */

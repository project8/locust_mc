/*
 * LMCButterworthLPFGenerator.cc
 *
 *  Created on: Sept 9, 2016
 *      Author: plslocum after nsoblath
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
            fDoGenerateFunc( &ButterworthLPFGenerator::DoGenerateTime )
    {
        fRequiredSignalState = Signal::kTime;
    }

    ButterworthLPFGenerator::~ButterworthLPFGenerator()
    {
    }

    bool ButterworthLPFGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;
        return true;
    }

    void ButterworthLPFGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    bool ButterworthLPFGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool ButterworthLPFGenerator::DoGenerateTime( Signal* aSignal ) const
    {
    	// 8th order Butterworth filter with wc = 70.e6 Hz, 100X attenuation at 95 MHz, fs=200 MHz.

    	double* ai = new double[20];
    	double* bi = new double[20];
    	int N = 8;  // 8th order digital filter.
    	int M = 8;

// raw coefficients from numerator of G(z) in Mathematica, starting with Butterworth polynomials.
    	ai[0] = 12.675;
    	ai[1] = 101.4;
    	ai[2] = 354.9;
    	ai[3] = 709.8;
    	ai[4] = 887.25;
    	ai[5] = 709.8;
    	ai[6] = 354.9;
    	ai[7] = 101.4;
    	ai[8] = 12.675;

// raw coefficients from denominator of G(z) in Mathematica, starting with Butterworth polynomials.
    	bi[0] = 160.698;
    	bi[1] = -511.67;
    	bi[2] = -832.628;
    	bi[3] = -837.901;
    	bi[4] = -561.456;
    	bi[5] = -252.64;
    	bi[6] = -73.9998;
    	bi[7] = -12.8136;
    	bi[8] = -1.;


double norm = bi[0];
for (int i=0; i<N+1; i++)
    {
	ai[i] /= norm;
	bi[i] /= norm;
    }

for (int i=0; i<N+1; i++)
  printf("ai[%d] is %f and bi[%d] is %f\n", i, ai[i], i, bi[i]);

    double FIR_component = 0.;
    double IIR_component = 0.;
    double *FilteredSignal = new double[aSignal->TimeSize()]; // temporary signal.
    for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
    	FilteredSignal[index] = 0.;


    for( unsigned index = N; index < aSignal->TimeSize(); ++index )
      {

      FIR_component = 0.;
      IIR_component = 0.;

      for (int i=0; i<N+1; i++)
	    FIR_component += ai[i]*aSignal->SignalTime()[index-i];
	  for (int i=1; i<M+1; i++)
	    IIR_component += bi[i]*FilteredSignal[index-i];

	  FilteredSignal[index] = FIR_component + IIR_component;
//      printf("FilteredSignal[%d] is %g\n", index, FilteredSignal[index]); getchar();
      }

    for( unsigned index = N+1; index < aSignal->TimeSize(); ++index )
      {
//      printf("FilteredSignal[%d] is %g\n", index, FilteredSignal[index]); getchar();
      aSignal->SignalTime()[index] = FilteredSignal[index];
      }


    delete FilteredSignal;
    return true;
    }

    bool ButterworthLPFGenerator::DoGenerateFreq( Signal* aSignal ) const
    {

        return true;
    }

} /* namespace locust */

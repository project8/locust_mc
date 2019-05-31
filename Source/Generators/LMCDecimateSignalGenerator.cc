/*
 * LMCDecimateSignalGenerator.cc
 *
 *  Created on: Sept 9, 2016
 *      Author: plslocum after nsoblath
 */

#include "LMCDecimateSignalGenerator.hh"

#include "logger.hh"


using std::string;

namespace locust
{
  LOGGER( lmclog, "DecimateSignalGenerator" );

  MT_REGISTER_GENERATOR(DecimateSignalGenerator, "decimate-signal");

  DecimateSignalGenerator::DecimateSignalGenerator( const std::string& aName ) :
    Generator( aName ),
    fDoGenerateFunc( &DecimateSignalGenerator::DoGenerateTime )
  {
    fRequiredSignalState = Signal::kTime;
  }

  DecimateSignalGenerator::~DecimateSignalGenerator()
  {
  }

  bool DecimateSignalGenerator::Configure( const scarab::param_node* aParam )
  {
    if( aParam == NULL) return true;
    return true;
  }

  void DecimateSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
  {
    aVisitor->Visit( this );
    return;
  }


  bool DecimateSignalGenerator::DoGenerate( Signal* aSignal )
  {
    return (this->*fDoGenerateFunc)( aSignal );
  }

  bool DecimateSignalGenerator::DoGenerateTime( Signal* aSignal ) 
  {
    const unsigned nchannels = fNChannels;

    // create text file for VI and VQ for testing. Assumes 1 channel.
    std::ofstream voltagefile;
    voltagefile.open("voltagefile.txt");
	
    // Decimate Fs -> Fs/10
    for (int ch=0; ch<nchannels; ch++)
      {
	for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
	  {
	    if (index % aSignal->DecimationFactor() == 0)
	      {
		aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index/aSignal->DecimationFactor()][0] =
		  aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][0];
		aSignal->SignalTimeComplex()[ch*aSignal->TimeSize() + index/aSignal->DecimationFactor()][1] =
		  aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + index][1];
	      }
	    // Save pre-decimated voltages for testing
	    voltagefile << index;
	    voltagefile << "\n";
	    voltagefile << aSignal->LongSignalTimeComplex()[index][0];
	    voltagefile << "\n";
	    voltagefile << aSignal->LongSignalTimeComplex()[index][1];
	    voltagefile << "\n";
	  }
      }
    voltagefile.close();

    return true;
  }

} /* namespace locust */

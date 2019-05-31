/*
 * LMCTestFIRFilter.cc
 *
 *  Created on: May 3 2019
 *      Author: plslocum
 */

#include <cmath>
#include "LMCTestFIRFilterGenerator.hh"
#include "LMCDigitizer.hh"
#include "logger.hh"
#include "LMCConst.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "TestFIRFilterGenerator" );

    MT_REGISTER_GENERATOR(TestFIRFilterGenerator, "test-firfilter");

    TestFIRFilterGenerator::TestFIRFilterGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &TestFIRFilterGenerator::DoGenerateTime ),
        gfilter_filename("blank.txt"),
        fFilter_resolution( 0. ),
        fLO_frequency( 20.05e9 ),
        fRF_frequency( 20.1e9 ),
        fAmplitude( 5.e-8 ),
        EFieldBuffer( 1 ),
        EPhaseBuffer( 1 ),
        EAmplitudeBuffer( 1 ),
        EFrequencyBuffer( 1 ),
        LOPhaseBuffer( 1 ),
        IndexBuffer( 1 ),
        PatchFIRBuffer( 1 ),
        fFieldBufferSize( 50 ),
        fFieldBufferMargin( 10 ),
		fNPatches( 1 )

    {
        fRequiredSignalState = Signal::kTime;
    }

    TestFIRFilterGenerator::~TestFIRFilterGenerator()
    {
    }


    bool TestFIRFilterGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        if( aParam->has( "filter-filename" ) )
        {
            gfilter_filename = aParam->get_value< std::string >( "filter-filename" );
        }

        if( aParam->has( "filter-resolution" ) )
        {
            fFilter_resolution = aParam->get_value< double >( "filter-resolution" );
        }


        if( aParam->has( "rf-frequency" ) )
        {
        SetRFFrequency( aParam->get_value< double >( "rf-frequency", fRF_frequency ) );
        }

        if( aParam->has( "lo-frequency" ) )
        {
        SetLOFrequency( aParam->get_value< double >( "lo-frequency", fLO_frequency ) );
        }

        if( aParam->has( "buffer-size" ) )
        {
        SetBufferSize( aParam->get_value< double >( "buffer-size", fFieldBufferSize ) );
        }

        if( aParam->has( "buffer-margin" ) )
        {
        SetBufferMargin( aParam->get_value< double >( "buffer-margin", fFieldBufferMargin ) );
        }

        if( aParam->has( "amplitude" ) )
        {
        SetAmplitude( aParam->get_value< double >( "amplitude", fAmplitude ) );
        }


        if( aParam->has( "domain" ) )
        {
            string domain = aParam->get_value( "domain" );
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                LDEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                LERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
                return false;
            }
        }


        return true;



    }


    void TestFIRFilterGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double TestFIRFilterGenerator::GetRFFrequency() const
    {
        return fRF_frequency;
    }

    void TestFIRFilterGenerator::SetRFFrequency( double aFrequency )
    {
        fRF_frequency = aFrequency;
        return;
    }

    double TestFIRFilterGenerator::GetLOFrequency() const
    {
        return fLO_frequency;
    }

    void TestFIRFilterGenerator::SetLOFrequency( double aFrequency )
    {
        fLO_frequency = aFrequency;
        return;
    }

    double TestFIRFilterGenerator::GetAmplitude() const
    {
        return fAmplitude;
    }

    void TestFIRFilterGenerator::SetAmplitude( double aAmplitude )
    {
        fAmplitude = aAmplitude;
        return;
    }

    double TestFIRFilterGenerator::GetBufferSize() const
    {
        return fFieldBufferSize;
    }

    void TestFIRFilterGenerator::SetBufferSize( double aBufferSize )
    {
        fFieldBufferSize = aBufferSize;
        return;
    }

    double TestFIRFilterGenerator::GetBufferMargin() const
    {
        return fFieldBufferMargin;
    }

    void TestFIRFilterGenerator::SetBufferMargin( double aBufferMargin )
    {
        fFieldBufferMargin = aBufferMargin;
        return;
    }



    Signal::State TestFIRFilterGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void TestFIRFilterGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &TestFIRFilterGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &TestFIRFilterGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    double* TestFIRFilterGenerator::GetFIRFilter(int nskips)
    {

    FILE *fp;
    double *filterarray = new double[1000];
    double filter;
    double index;
    fp = fopen(gfilter_filename.c_str(),"r");
    int count = 0;


    for (int i=0; i<1000; i++)
      filterarray[i] = -99.;



      while (!feof(fp))
        {
        fscanf(fp, "%lf %lf\n", &index, &filter);
        if (count%nskips==0) filterarray[count/nskips] = filter;
    //    printf("filter %d is %g\n", count, filterarray[count]);
        count += 1;
        }

    fclose(fp);
    return filterarray;

    }

    int TestFIRFilterGenerator::GetNFilterBins(double* filterarray)
    {
    int nbins = 0;
    for (int i=0; i<1000; i++)
      {
      if (filterarray[i]>0.) nbins += 1;
      }    
    return nbins;
    }


    double TestFIRFilterGenerator::GetFIRSample(double* filterarray, int nfilterbins, double dtfilter, unsigned channel, unsigned patch, double AcquisitionRate)
    {

    double fieldfrequency = EFrequencyBuffer[channel*fNPatches+patch].front();
    double convolution = 0.;
    double HilbertMag = 0.;
    double HilbertPhase = 0.;
    double HilbertMean = 0.;

    if (fabs(EFieldBuffer[channel].front()) > 0.)  // field arrived yet?
    {
        HilbertTransform aHilbertTransform;
        double* HilbertMagPhaseMean = aHilbertTransform.GetMagPhaseMean(EFieldBuffer[channel*fNPatches+patch], EFrequencyBuffer[channel*fNPatches+patch], fFieldBufferMargin, AcquisitionRate);
        HilbertMag = HilbertMagPhaseMean[0];
        HilbertPhase = HilbertMagPhaseMean[1];
        HilbertMean = HilbertMagPhaseMean[2];
//        HilbertMag = fAmplitude;  // cross check with analytic inputs.
//        HilbertPhase = EPhaseBuffer[channel].front(); // cross check with analytic inputs.


   	for (int i=0; i<nfilterbins; i++)  // populate filter with field.
    	  {
    	  HilbertPhase += 2.*3.1415926*fieldfrequency*dtfilter;
    	  PatchFIRBuffer[channel*fNPatches+patch].push_back(HilbertMag*cos(HilbertPhase));
    	  PatchFIRBuffer[channel*fNPatches+patch].pop_front();
    	  }

    	for (int j=0; j<nfilterbins; j++)  // sum products in filter.
    	  {
    	  convolution += PatchFIRBuffer[channel*fNPatches+patch].at(j)*filterarray[j];
    	  }

    }

//    return sqrt(50.)*HilbertMag*cos(HilbertPhase);  // debug
    return convolution;
    }


    bool TestFIRFilterGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool TestFIRFilterGenerator::DoGenerateTime( Signal* aSignal )
    {


        const unsigned nchannels = fNChannels;
        const unsigned npatches = fNPatches;

        double LO_phase = 0.;
        double field_phase = 0.;
        double VoltageSample = 0.;
        double* filterarray = GetFIRFilter(1);
        unsigned nfilterbins = GetNFilterBins(filterarray);
        unsigned nfieldbufferbins = fFieldBufferSize;
        double dtfilter = fFilter_resolution;

        InitializeBuffers(nfilterbins, nfieldbufferbins);

        for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
          {
          for (unsigned ch = 0; ch < nchannels; ++ch)
            {
        	for (unsigned patch = 0; patch < npatches; ++patch)
        	{

            LO_phase += 2.*LMCConst::Pi()*fLO_frequency/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
            field_phase += 2.*LMCConst::Pi()*fRF_frequency/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);

            FillBuffers(fAmplitude, field_phase, LO_phase, index, ch, patch);
            VoltageSample = GetFIRSample(filterarray, nfilterbins, dtfilter, ch, patch, fAcquisitionRate*aSignal->DecimationFactor());
//            VoltageSample = sqrt(50.)*fAmplitude*cos(field_phase);
//            printf("voltagesample is %g\n", VoltageSample); getchar();

// factor of 2 is needed for cosA*cosB = 1/2*(cos(A+B)+cos(A-B)); usually we leave out the 1/2 for e.g. sinusoidal RF.
// This allows for correct gain in Locust-Katydid analysis chain.
            aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + (unsigned)IndexBuffer[ch].front()][0] += 2.*VoltageSample*cos(LOPhaseBuffer[ch].front());  
            aSignal->LongSignalTimeComplex()[ch*aSignal->TimeSize()*aSignal->DecimationFactor() + (unsigned)IndexBuffer[ch].front()][1] += 2.*VoltageSample*cos(LMCConst::Pi()/2. + LOPhaseBuffer[ch].front());

            PopBuffers(ch, patch);

        	}  // npatches

            }  // nchannels
        }  // index

        return true;
    }

    void TestFIRFilterGenerator::FillBuffers(double FieldAmplitude, double FieldPhase, double LOPhase, unsigned index, unsigned channel, unsigned patch)
    {
    EFieldBuffer[channel*fNPatches+patch].push_back(fAmplitude*cos(FieldPhase));
    EPhaseBuffer[channel*fNPatches+patch].push_back(FieldPhase);
    EAmplitudeBuffer[channel*fNPatches+patch].push_back(fAmplitude);
    EFrequencyBuffer[channel*fNPatches+patch].push_back(fRF_frequency);
    LOPhaseBuffer[channel*fNPatches+patch].push_back(LOPhase);
    IndexBuffer[channel*fNPatches+patch].push_back(index);
    }


    void TestFIRFilterGenerator::PopBuffers(unsigned channel, unsigned patch)
    {
    EFieldBuffer[channel*fNPatches+patch].pop_front();
    EPhaseBuffer[channel*fNPatches+patch].pop_front();
    EAmplitudeBuffer[channel*fNPatches+patch].pop_front();
    EFrequencyBuffer[channel*fNPatches+patch].pop_front();
    LOPhaseBuffer[channel*fNPatches+patch].pop_front();
    IndexBuffer[channel*fNPatches+patch].pop_front();
    }


    void TestFIRFilterGenerator::InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize)
    {

    const unsigned nchannels = fNChannels;

    FieldBuffer aFieldBuffer;

    EFieldBuffer = aFieldBuffer.InitializeBuffer(nchannels, 1, fieldbuffersize);
    EPhaseBuffer = aFieldBuffer.InitializeBuffer(nchannels, 1, fieldbuffersize);
    EAmplitudeBuffer = aFieldBuffer.InitializeBuffer(nchannels, 1, fieldbuffersize);
    EFrequencyBuffer = aFieldBuffer.InitializeBuffer(nchannels, 1, fieldbuffersize);
    LOPhaseBuffer = aFieldBuffer.InitializeBuffer(nchannels, 1, fieldbuffersize);
    IndexBuffer = aFieldBuffer.InitializeBuffer(nchannels, 1, fieldbuffersize);

    PatchFIRBuffer = aFieldBuffer.InitializeBuffer(1, 1, filterbuffersize);
    }

    bool TestFIRFilterGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */

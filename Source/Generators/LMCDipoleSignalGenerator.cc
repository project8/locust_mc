/*
 * LMCDipoleSignalGenerator.cc
 *
 *  Created on: July 11 2019
 *      Author: Pranava Teja Surukuchi
 */

#include <cmath>
#include "LMCDipoleSignalGenerator.hh"
#include "LMCDigitizer.hh"
#include "logger.hh"
#include "LMCConst.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "DipoleSignalGenerator" );

    MT_REGISTER_GENERATOR(DipoleSignalGenerator, "antenna-dipole-generator");

    DipoleSignalGenerator::DipoleSignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fDoGenerateFunc( &DipoleSignalGenerator::DoGenerateTime ),
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
	ConvolutionTimeBuffer( 1 ),
        fFieldBufferSize( 50 ),
        fFieldBufferMargin( 25 ),
	fNPatches( 1 )

    {
        fRequiredSignalState = Signal::kTime;
    }

    DipoleSignalGenerator::~DipoleSignalGenerator()
    {
    }


    bool DipoleSignalGenerator::Configure( const scarab::param_node* aParam )
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


        if( aParam->has( "input-signal-frequency"))
        {
        SetRFFrequency( aParam->get_value< double >( "input-signal-frequency", fRF_frequency ) );
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

        if( aParam->has( "input-signal-amplitude" ) )
        {
        SetAmplitude( aParam->get_value< double >( "input-signal-amplitude", fAmplitude ) );
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

	if(!fAntennaSignalGenerator.Configure(aParam))
	{
		LERROR(lmclog,"Error Configuring antenna signal generator class");
	}	
        return true;
    }


    void DipoleSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    double DipoleSignalGenerator::GetRFFrequency() const
    {
        return fRF_frequency;
    }

    void DipoleSignalGenerator::SetRFFrequency( double aFrequency )
    {
        fRF_frequency = aFrequency;
        return;
    }

    double DipoleSignalGenerator::GetLOFrequency() const
    {
        return fLO_frequency;
    }

    void DipoleSignalGenerator::SetLOFrequency( double aFrequency )
    {
        fLO_frequency = aFrequency;
        return;
    }

    double DipoleSignalGenerator::GetAmplitude() const
    {
        return fAmplitude;
    }

    void DipoleSignalGenerator::SetAmplitude( double aAmplitude )
    {
        fAmplitude = aAmplitude;
        return;
    }

    double DipoleSignalGenerator::GetBufferSize() const
    {
        return fFieldBufferSize;
    }

    void DipoleSignalGenerator::SetBufferSize( double aBufferSize )
    {
        fFieldBufferSize = aBufferSize;
        return;
    }

    double DipoleSignalGenerator::GetBufferMargin() const
    {
        return fFieldBufferMargin;
    }

    void DipoleSignalGenerator::SetBufferMargin( double aBufferMargin )
    {
        fFieldBufferMargin = aBufferMargin;
        return;
    }



    Signal::State DipoleSignalGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void DipoleSignalGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &DipoleSignalGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &DipoleSignalGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    double* DipoleSignalGenerator::GetFIRFilter(int nskips)
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

    int DipoleSignalGenerator::GetNFilterBins(double* filterarray)
    {
    int nbins = 0;
    for (int i=0; i<1000; i++)
      {
      if (filterarray[i]>0.) nbins += 1;
      }    
    return nbins;
    }


    double DipoleSignalGenerator::GetFIRSample(double* filterarray, int nfilterbins, double dtfilter, unsigned channel, unsigned patch, double AcquisitionRate)
    {

    double fieldfrequency = EFrequencyBuffer[channel*fNPatches+patch].front();
    double HilbertMag = 0.;
    double HilbertPhase = 0.;
    double convolution = 0.0;

    if (fabs(EFieldBuffer[channel*fNPatches+patch].front()) > 0.)  // field arrived yet?
    {
    	HilbertTransform aHilbertTransform;

    	double* HilbertMagPhaseMean = new double[3];
        HilbertMagPhaseMean = aHilbertTransform.GetMagPhaseMean(EFieldBuffer[channel*fNPatches+patch], EFrequencyBuffer[channel*fNPatches+patch], fFieldBufferMargin, AcquisitionRate);
        HilbertMag = HilbertMagPhaseMean[0];
        HilbertPhase = HilbertMagPhaseMean[1];
        delete[] HilbertMagPhaseMean;

   	for (int i=0; i < nfilterbins - ConvolutionTimeBuffer[channel*fNPatches+patch].front(); i++)  // populate filter with field.
      {
    	  HilbertPhase += 2.*LMCConst::Pi()*fieldfrequency*dtfilter;
    	  PatchFIRBuffer[channel*fNPatches+patch].push_back(HilbertMag*cos(HilbertPhase));
    	  PatchFIRBuffer[channel*fNPatches+patch].pop_front();
      }

    for (int j=0; j<nfilterbins; j++)  // sum products in filter.
      {
    	  convolution += filterarray[j]*PatchFIRBuffer[channel*fNPatches+patch].at(j);
      }

    PatchFIRBuffer[channel*fNPatches+patch].shrink_to_fit();  // memory deallocation.
    return convolution;
    }
    else return 0.;

//    return sqrt(50.)*HilbertMag*cos(HilbertPhase);  // debug
//    printf("returning convolution %d %d\n", channel, patch);
    }

    bool DipoleSignalGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }


    bool DipoleSignalGenerator::DoGenerateTime( Signal* aSignal )
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
        unsigned dtauConvolutionTime = 0;

	fAntennaSignalGenerator.InitializeTransmitter();
        InitializeBuffers(nfilterbins, nfieldbufferbins);

	field_phase=fAntennaSignalGenerator.GetInitialPhaseDelay();	
        for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
        {
		double fieldValue=fAntennaSignalGenerator.GenerateSignal(aSignal,fAcquisitionRate);
          	LO_phase += 2.*LMCConst::Pi()*fLO_frequency/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
          	field_phase += 2.*LMCConst::Pi()*fRF_frequency/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6) ;
          	for (unsigned ch = 0; ch < nchannels; ++ch)
            	{
        		for (unsigned patch = 0; patch < npatches; ++patch)
        		{
            			if (index > 0) dtauConvolutionTime = 0;
            			else dtauConvolutionTime = nfilterbins/2;
            			FillBuffers(aSignal, fieldValue, field_phase, LO_phase, index, ch, patch, dtauConvolutionTime);
		            	VoltageSample = GetFIRSample(filterarray, nfilterbins, dtfilter, ch, patch, fAcquisitionRate*aSignal->DecimationFactor());
// factor of 2 is needed for cosA*cosB = 1/2*(cos(A+B)+cos(A-B)); usually we leave out the 1/2 for e.g. sinusoidal RF.
            			aSignal->LongSignalTimeComplex()[IndexBuffer[ch*fNPatches+patch].front()][0] += 2.*VoltageSample*cos(LOPhaseBuffer[ch*fNPatches+patch].front());
            			aSignal->LongSignalTimeComplex()[IndexBuffer[ch*fNPatches+patch].front()][1] += 2.*VoltageSample*(-sin(LOPhaseBuffer[ch*fNPatches+patch].front()));
				PopBuffers(ch, patch);
        		}  // patch
            	}  // channel
		//std::cout<< "index: "<< index <<std::endl;
        }  // index
        CleanupBuffers();
        return true;
    }

    void DipoleSignalGenerator::FillBuffers(Signal* aSignal, double FieldValue, double FieldPhase, double LOPhase, unsigned index, unsigned channel, unsigned patch, unsigned dtauConvolutionTime)
    {

    EFieldBuffer[channel*fNPatches+patch].push_back(FieldValue);
    EPhaseBuffer[channel*fNPatches+patch].push_back(FieldPhase);
    EAmplitudeBuffer[channel*fNPatches+patch].push_back(fAmplitude);
    EFrequencyBuffer[channel*fNPatches+patch].push_back(fRF_frequency);
    LOPhaseBuffer[channel*fNPatches+patch].push_back(LOPhase);
    IndexBuffer[channel*fNPatches+patch].push_back(channel*aSignal->TimeSize()*aSignal->DecimationFactor() + index);
    ConvolutionTimeBuffer[channel*fNPatches+patch].push_back(dtauConvolutionTime);

    }


    void DipoleSignalGenerator::PopBuffers(unsigned channel, unsigned patch)
    {
    	EFieldBuffer[channel*fNPatches+patch].pop_front();
    	EPhaseBuffer[channel*fNPatches+patch].pop_front();
    	EAmplitudeBuffer[channel*fNPatches+patch].pop_front();
    	EFrequencyBuffer[channel*fNPatches+patch].pop_front();
    	LOPhaseBuffer[channel*fNPatches+patch].pop_front();
    	IndexBuffer[channel*fNPatches+patch].pop_front();
    	ConvolutionTimeBuffer[channel*fNPatches+patch].pop_front();

    	EFieldBuffer[channel*fNPatches+patch].shrink_to_fit();
        EPhaseBuffer[channel*fNPatches+patch].shrink_to_fit();
        EAmplitudeBuffer[channel*fNPatches+patch].shrink_to_fit();
        EFrequencyBuffer[channel*fNPatches+patch].shrink_to_fit();
        LOPhaseBuffer[channel*fNPatches+patch].shrink_to_fit();
        IndexBuffer[channel*fNPatches+patch].shrink_to_fit();
        // PTS: Seg faults when shrink-to-fit used on ConvolutionTimeBuffer, removed for now. Need to revisit what the problem is
	//ConvolutionTimeBuffer[channel+fNPatches+patch].shrink_to_fit();

    }




    void DipoleSignalGenerator::CleanupBuffers()
    {
    FieldBuffer aFieldBuffer;
    EFieldBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    EPhaseBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    EAmplitudeBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    EFrequencyBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    LOPhaseBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    IndexBuffer = aFieldBuffer.CleanupBuffer(IndexBuffer);
    ConvolutionTimeBuffer = aFieldBuffer.CleanupBuffer(ConvolutionTimeBuffer);

    }

    void DipoleSignalGenerator::InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize)
    {

    FieldBuffer aFieldBuffer;

    EFieldBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatches, fieldbuffersize);
    EPhaseBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatches, fieldbuffersize);
    EAmplitudeBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatches, fieldbuffersize);
    EFrequencyBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatches, fieldbuffersize);
    LOPhaseBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatches, fieldbuffersize);
    IndexBuffer = aFieldBuffer.InitializeUnsignedBuffer(fNChannels, fNPatches, fieldbuffersize);
    ConvolutionTimeBuffer = aFieldBuffer.InitializeUnsignedBuffer(fNChannels, fNPatches, fieldbuffersize);

    PatchFIRBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatches, filterbuffersize);


    }


    bool DipoleSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }

} /* namespace locust */

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
    fArrayRadius( 0. ),
    fNPatchesPerStrip( 1 ),
    fPatchSpacing( 0. ),
    fPowerCombiner( 0 ),
    fTextFileWriting( 0 ),
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
    fFieldBufferMargin( 25 )
    
    {
        fRequiredSignalState = Signal::kTime;
    }
    
    DipoleSignalGenerator::~DipoleSignalGenerator()
    {
    }
    
    
    bool DipoleSignalGenerator::Configure( const scarab::param_node& aParam )
    {
        if( aParam.has( "filter-filename" ) )
        {
            gfilter_filename = aParam["filter-filename"]().as_string();
        }
        
        if( aParam.has( "filter-resolution" ) )
        {
            fFilter_resolution = aParam["filter-resolution"]().as_double();
        }
        
        if( aParam.has( "array-radius" ) )
        {
            fArrayRadius = aParam["array-radius"]().as_double();
        }
        
        if( aParam.has( "input-signal-frequency"))
        {
            SetRFFrequency( aParam.get_value< double >( "input-signal-frequency", fRF_frequency ) );
        }
        
        if( aParam.has( "lo-frequency" ) )
        {
            SetLOFrequency( aParam.get_value< double >( "lo-frequency", fLO_frequency ) );
        }
        
        if( aParam.has( "buffer-size" ) )
        {
            SetBufferSize( aParam.get_value< double >( "buffer-size", fFieldBufferSize ) );
        }
        
        if( aParam.has( "buffer-margin" ) )
        {
            SetBufferMargin( aParam.get_value< double >( "buffer-margin", fFieldBufferMargin ) );
        }
        
        if( aParam.has( "input-signal-amplitude" ) )
        {
            SetAmplitude( aParam.get_value< double >( "input-signal-amplitude", fAmplitude ) );
        }
        
        if( aParam.has( "npatches-per-strip" ) )
        {
            fNPatchesPerStrip = aParam["npatches-per-strip"]().as_int();
        }
        
        if( aParam.has( "patch-spacing" ) )
        {
            fPatchSpacing = aParam["patch-spacing"]().as_double();
        }
        if( aParam.has( "xml-filename" ) )
        {
            gxml_filename = aParam["xml-filename"]().as_string();
        }
        if( aParam.has( "text-filewriting" ) )
        {
            fTextFileWriting = aParam["text-filewriting"]().as_bool();
        }
        if( aParam.has( "feed" ) )
        {
            if (aParam["feed"]().as_string() == "corporate")
                fPowerCombiner = 0;  // default
            else if (aParam["feed"]().as_string() == "series")
                fPowerCombiner = 1;
            else if (aParam["feed"]().as_string() == "one-quarter")
                fPowerCombiner = 2;
            else if (aParam["feed"]().as_string() == "seven-eighths")
                fPowerCombiner = 3;
            else if (aParam["feed"]().as_string() == "nine-sixteenths")
                fPowerCombiner = 4;
            else
                fPowerCombiner = 0;  // default
        }
        
        if( aParam.has( "domain" ) )
        {
            string domain = aParam["domain"]().as_string();
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
        
        if(!fAntennaSignalTransmitter.Configure(aParam))
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
    
    //PTS: This part(function) should be unified with the other parts where a filter is extracted
    double DipoleSignalGenerator::GetAOIFactor(LMCThreeVector TrasmittingPatchNormal,LMCThreeVector ReceivingPatchNormal)
    {
        double AOIFactor = fabs(TrasmittingPatchNormal.Unit().Dot(ReceivingPatchNormal.Unit()));
        return AOIFactor;
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
    
    //PTS: This function should be unified with the other FIR Filter-using generators where a filter is extracted
    double DipoleSignalGenerator::GetFIRSample(double* filterarray, int nfilterbins, double dtfilter, unsigned channel, unsigned patch,double fieldPhase, double AcquisitionRate)
    {
        
        double fieldfrequency = EFrequencyBuffer[channel*fNPatchesPerStrip+patch].front();
        double HilbertMag = 0.;
        double HilbertPhase = 0.;
        double convolution = 0.0;
        
        if (fabs(EFieldBuffer[channel*fNPatchesPerStrip+patch].front()) > 0.)  // field arrived yet?
        {
            HilbertTransform aHilbertTransform;
            
            double* HilbertMagPhaseMean = new double[3];
            HilbertMagPhaseMean = aHilbertTransform.GetMagPhaseMean(EFieldBuffer[channel*fNPatchesPerStrip+patch], EFrequencyBuffer[channel*fNPatchesPerStrip+patch], fFieldBufferMargin, AcquisitionRate);
            HilbertMag = HilbertMagPhaseMean[0];
            HilbertPhase = HilbertMagPhaseMean[1];
            delete[] HilbertMagPhaseMean;
            
            HilbertPhase += fieldPhase;
            for (int i=0; i < nfilterbins - ConvolutionTimeBuffer[channel*fNPatchesPerStrip+patch].front(); i++)  // populate filter with field.
            {
                HilbertPhase += 2.*LMCConst::Pi()*fieldfrequency*dtfilter;
                PatchFIRBuffer[channel*fNPatchesPerStrip+patch].push_back(HilbertMag*cos(HilbertPhase));
                PatchFIRBuffer[channel*fNPatchesPerStrip+patch].pop_front();
            }
            
            for (int j=0; j<nfilterbins; j++)  // sum products in filter.
            {
                convolution += filterarray[j]*PatchFIRBuffer[channel*fNPatchesPerStrip+patch].at(j);
            }
            
            PatchFIRBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();  // memory deallocation.
            return convolution;
        }
        else return 0.;
        
    }
    
    // Should be moved to somewhere more general
    double DipoleSignalGenerator::RotateZ(int component, double angle, double x, double y)
    {
        double newcomponent = 0.;
        if (component==0)
        {
            newcomponent = x*cos(angle) - y*sin(angle);
        }
        else if (component==1)
        {
            newcomponent = x*sin(angle) + y*cos(angle);
        }
        return newcomponent;
    }
    
    void DipoleSignalGenerator::AddOneFIRVoltageToStripSum(Signal* aSignal, double VoltageSample, double phi_LO, unsigned channelIndex, unsigned patchIndex)
    {
        
        PowerCombiner aPowerCombiner;
        
        if (fPowerCombiner == 0 ) //corporate feed, for testing
        {
            VoltageSample *= aPowerCombiner.GetCorporateVoltageDamping();
        }
        
        if (fPowerCombiner == 3) // seven-eighths power combining, center fed strip
        {
            VoltageSample *= aPowerCombiner.GetSevenEighthsVoltageDamping(fNPatchesPerStrip, patchIndex);
        }
        
        aSignal->LongSignalTimeComplex()[IndexBuffer[channelIndex*fNPatchesPerStrip+patchIndex].front()][0] += 2.*VoltageSample * cos(phi_LO);
        aSignal->LongSignalTimeComplex()[IndexBuffer[channelIndex*fNPatchesPerStrip+patchIndex].front()][1] += 2.*VoltageSample * sin(phi_LO);
        
    }
    
    void  DipoleSignalGenerator::InitializePatchArray()
    {
        
        const unsigned nChannels = fNChannels;
        const int nReceivers = fNPatchesPerStrip;
        const double patchSpacingZ = fPatchSpacing;
        const double patchRadius = fArrayRadius;
        double zPosition;
        double theta;
        const double dThetaArray = 2. * LMCConst::Pi() / nChannels; //Divide the circle into nChannels
        const double dRotateVoltages = 0.;  // set to zero to not rotate patch polarities.
        
        PatchAntenna modelPatch;
        
        allChannels.resize(nChannels);
        
        for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        {
            theta = channelIndex * dThetaArray;
            
            for(int receiverIndex = 0; receiverIndex < nReceivers; ++receiverIndex)
            {
                zPosition =  (receiverIndex - (nReceivers - 1.) /2.) * patchSpacingZ;
                
                modelPatch.SetCenterPosition({patchRadius * cos(theta) , patchRadius * sin(theta) , zPosition });
                modelPatch.SetPolarizationDirection({RotateZ(0, dRotateVoltages*channelIndex, sin(theta), -cos(theta)), RotateZ(1, dRotateVoltages*channelIndex, sin(theta), -cos(theta)), 0.});
                
                modelPatch.SetNormalDirection({-cos(theta), -sin(theta), 0.}); //Say normals point inwards
                allChannels[channelIndex].AddReceiver(modelPatch);
            }
        }
    }
    
    bool DipoleSignalGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }
    
    
    bool DipoleSignalGenerator::DoGenerateTime( Signal* aSignal )
    {
        
        InitializePatchArray();
        const unsigned nchannels = fNChannels;
        const unsigned npatches = fNPatchesPerStrip;
        
        double LO_phase = 0.;
        double field_phase = 0.;
        double VoltageSample = 0.;
        double* filterarray = GetFIRFilter(1);
        unsigned nfilterbins = GetNFilterBins(filterarray);
        unsigned nfieldbufferbins = fFieldBufferSize;
        double dtfilter = fFilter_resolution;
        unsigned dtauConvolutionTime = 0;
        
        if(!fAntennaSignalTransmitter.InitializeTransmitter())
        {
            exit(-1);
        }
        
        InitializeBuffers(nfilterbins, nfieldbufferbins);
        
        double initialPhaseDelay=fAntennaSignalTransmitter.GetInitialPhaseDelay();
        for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
        {
            double fieldValue=fAntennaSignalTransmitter.GenerateSignal(aSignal,fAcquisitionRate);
            double antennaPositionX=fAntennaSignalTransmitter.GetAntennaPosition().GetX();
            double antennaPositionY=fAntennaSignalTransmitter.GetAntennaPosition().GetY();
            double antennaPositionZ=fAntennaSignalTransmitter.GetAntennaPosition().GetZ();
            LO_phase += 2.*LMCConst::Pi()*fLO_frequency/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
            for (unsigned ch = 0; ch < nchannels; ++ch)
            {
                for (unsigned patch = 0; patch < npatches; ++patch)
                {
                    PatchAntenna *currentPatch;
                    currentPatch = &allChannels[ch][patch];
                    double relativePatchPosX=currentPatch->GetPosition().GetX() - antennaPositionX;
                    double relativePatchPosY=currentPatch->GetPosition().GetY() - antennaPositionY;
                    double relativePatchPosZ=currentPatch->GetPosition().GetZ() - antennaPositionZ;
                    double patchAntennaDistance = sqrt(relativePatchPosX*relativePatchPosX+relativePatchPosY*relativePatchPosY+relativePatchPosZ*relativePatchPosZ);
                    fieldValue=fieldValue*GetAOIFactor(fAntennaSignalTransmitter.GetAntennaPosition()-currentPatch->GetNormalDirection(),currentPatch->GetNormalDirection());
                    double field_phase=initialPhaseDelay+2.*LMCConst::Pi()*(patchAntennaDistance/LMCConst::C())*fRF_frequency;
                    if (index > 0) dtauConvolutionTime = 0;
                    else dtauConvolutionTime = nfilterbins/2;
                    FillBuffers(aSignal, fieldValue, field_phase, LO_phase, index, ch, patch, dtauConvolutionTime);
                    VoltageSample = GetFIRSample(filterarray, nfilterbins, dtfilter, ch, patch, field_phase,fAcquisitionRate*aSignal->DecimationFactor());
                    VoltageSample = VoltageSample/patchAntennaDistance;
                    AddOneFIRVoltageToStripSum(aSignal, VoltageSample, LO_phase, ch, patch);
                    // factor of 2 is needed for cosA*cosB = 1/2*(cos(A+B)+cos(A-B)); usually we leave out the 1/2 for e.g. sinusoidal RF.
                    //aSignal->LongSignalTimeComplex()[IndexBuffer[ch*fNPatchesPerStrip+patch].front()][0] += 2.*VoltageSample*cos(LOPhaseBuffer[ch*fNPatchesPerStrip+patch].front());
                    //aSignal->LongSignalTimeComplex()[IndexBuffer[ch*fNPatchesPerStrip+patch].front()][1] += 2.*VoltageSample*(-sin(LOPhaseBuffer[ch*fNPatchesPerStrip+patch].front()));
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
        
        EFieldBuffer[channel*fNPatchesPerStrip+patch].push_back(FieldValue);
        EPhaseBuffer[channel*fNPatchesPerStrip+patch].push_back(FieldPhase);
        EAmplitudeBuffer[channel*fNPatchesPerStrip+patch].push_back(fAmplitude);
        EFrequencyBuffer[channel*fNPatchesPerStrip+patch].push_back(fRF_frequency);
        LOPhaseBuffer[channel*fNPatchesPerStrip+patch].push_back(LOPhase);
        IndexBuffer[channel*fNPatchesPerStrip+patch].push_back(channel*aSignal->TimeSize()*aSignal->DecimationFactor() + index);
        ConvolutionTimeBuffer[channel*fNPatchesPerStrip+patch].push_back(dtauConvolutionTime);
        
    }
    
    
    void DipoleSignalGenerator::PopBuffers(unsigned channel, unsigned patch)
    {
        EFieldBuffer[channel*fNPatchesPerStrip+patch].pop_front();
        EPhaseBuffer[channel*fNPatchesPerStrip+patch].pop_front();
        EAmplitudeBuffer[channel*fNPatchesPerStrip+patch].pop_front();
        EFrequencyBuffer[channel*fNPatchesPerStrip+patch].pop_front();
        LOPhaseBuffer[channel*fNPatchesPerStrip+patch].pop_front();
        IndexBuffer[channel*fNPatchesPerStrip+patch].pop_front();
        ConvolutionTimeBuffer[channel*fNPatchesPerStrip+patch].pop_front();
        
        EFieldBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        EPhaseBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        EAmplitudeBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        EFrequencyBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        LOPhaseBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        IndexBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        // PTS: Seg faults when shrink-to-fit used on ConvolutionTimeBuffer, removed for now. Need to revisit what the problem is
        //ConvolutionTimeBuffer[channel+fNPatchesPerStrip+patch].shrink_to_fit();
        
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
        
        EFieldBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
        EPhaseBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
        EAmplitudeBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
        EFrequencyBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
        LOPhaseBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
        IndexBuffer = aFieldBuffer.InitializeUnsignedBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
        ConvolutionTimeBuffer = aFieldBuffer.InitializeUnsignedBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
        
        PatchFIRBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, filterbuffersize);
        
        
    }
    
    
    bool DipoleSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }
    
} /* namespace locust */

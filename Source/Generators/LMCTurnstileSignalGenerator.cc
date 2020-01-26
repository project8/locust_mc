/*
 * LMCTurnstileSignalGenerator.cc
 *
 *  Created on: Oct 2019
 *      Author: Pranava Teja Surukuchi
 */

#include <cmath>
#include "LMCTurnstileSignalGenerator.hh"
#include "LMCDigitizer.hh"
#include "logger.hh"
#include "LMCConst.hh"


using std::string;

namespace locust
{
    LOGGER( lmclog, "TurnstileSignalGenerator" );
    
    MT_REGISTER_GENERATOR(TurnstileSignalGenerator, "turnstile-signal-generator");
    
    TurnstileSignalGenerator::TurnstileSignalGenerator( const std::string& aName ) :
    Generator( aName ),
    fDoGenerateFunc( &TurnstileSignalGenerator::DoGenerateTime ),
    fLO_frequency( 20.05e9 ),
    fRF_frequency( 20.1e9 ),
    fArrayRadius( 0. ),
    fNElementsPerStrip( 1 ),
    fElementSpacing( 0. ),
    fTextFileWriting( 0 ),
    fAmplitude( 5.e-8 ),
    EFieldBuffer( 1 ),
    EPhaseBuffer( 1 ),
    EAmplitudeBuffer( 1 ),
    EFrequencyBuffer( 1 ),
    LOPhaseBuffer( 1 ),
    IndexBuffer( 1 ),
    PatchFIRBuffer( 1 ),
    fFieldBufferSize( 50 ),
	fSwapFrequency( 1000 )

    
    {
        fRequiredSignalState = Signal::kTime;
    }
    
    TurnstileSignalGenerator::~TurnstileSignalGenerator()
    {
    }
    
    
    bool TurnstileSignalGenerator::Configure( const scarab::param_node& aParam )
    {
        if(!fReceiverHandler.Configure(aParam))
        {
            LERROR(lmclog,"Error configuring receiver FIRHandler class");
        }

        if(!fPowerCombiner.Configure(aParam))
        {
            LERROR(lmclog,"Error configuring receiver FIRHandler class");
        }

        if( aParam.has( "buffer-size" ) )
        {
            SetBufferSize( aParam.get_value< double >( "buffer-size", fFieldBufferSize ) );
        	fHilbertTransform.SetBufferSize(aParam["buffer-size"]().as_int());
        }

    	if(!fHilbertTransform.Configure(aParam))
    	{
    		LERROR(lmclog,"Error configuring receiver HilbertTransform class");
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
        
        if( aParam.has( "input-signal-amplitude" ) )
        {
            SetAmplitude( aParam.get_value< double >( "input-signal-amplitude", fAmplitude ) );
        }
        
        if( aParam.has( "nelements-per-strip" ) )
        {
            fNElementsPerStrip = aParam["nelements-per-strip"]().as_int();
        }
        
        if( aParam.has( "element-spacing" ) )
        {
            fElementSpacing = aParam["element-spacing"]().as_double();
        }
        if( aParam.has( "xml-filename" ) )
        {
            gxml_filename = aParam["xml-filename"]().as_string();
        }
        if( aParam.has( "text-filewriting" ) )
        {
            fTextFileWriting = aParam["text-filewriting"]().as_bool();
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
    
    
    void TurnstileSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }
    
    double TurnstileSignalGenerator::GetRFFrequency() const
    {
        return fRF_frequency;
    }
    
    void TurnstileSignalGenerator::SetRFFrequency( double aFrequency )
    {
        fRF_frequency = aFrequency;
        return;
    }
    
    double TurnstileSignalGenerator::GetLOFrequency() const
    {
        return fLO_frequency;
    }
    
    void TurnstileSignalGenerator::SetLOFrequency( double aFrequency )
    {
        fLO_frequency = aFrequency;
        return;
    }
    
    double TurnstileSignalGenerator::GetAmplitude() const
    {
        return fAmplitude;
    }
    
    void TurnstileSignalGenerator::SetAmplitude( double aAmplitude )
    {
        fAmplitude = aAmplitude;
        return;
    }
    
    double TurnstileSignalGenerator::GetBufferSize() const
    {
        return fFieldBufferSize;
    }
    
    void TurnstileSignalGenerator::SetBufferSize( double aBufferSize )
    {
        fFieldBufferSize = aBufferSize;
        return;
    }
    
    
    Signal::State TurnstileSignalGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }
    
    void TurnstileSignalGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &TurnstileSignalGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &TurnstileSignalGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }
    
    //PTS: This part(function) should be unified with the other parts where a filter is extracted
    double TurnstileSignalGenerator::GetAOIFactor(LMCThreeVector TrasmitterReceiverNormal,LMCThreeVector ReceivingPatchNormal)
    {
	if(TrasmitterReceiverNormal==LMCThreeVector::sZero)
	{
	    return 1.0; //Unusual case where source is same as the receiver
	}
	double AOIFactor = fabs(TrasmitterReceiverNormal.Unit().Dot(ReceivingPatchNormal.Unit()));
        return AOIFactor;
    }
    
    double TurnstileSignalGenerator::GetVoltageFromField(unsigned channel, unsigned patch,double fieldPhase)
    {
        
        double fieldfrequency = EFrequencyBuffer[channel*fNElementsPerStrip+patch].front();
        double HilbertMag = 0.;
        double HilbertPhase = 0.;
        
        if (fabs(EFieldBuffer[channel*fNElementsPerStrip+patch].front()) > 0.)  // field arrived yet?
        {
            
            double* HilbertMagPhaseMean = new double[3];
            HilbertMagPhaseMean = fHilbertTransform.GetMagPhaseMean(EFieldBuffer[channel*fNElementsPerStrip+patch], EFrequencyBuffer[channel*fNElementsPerStrip+patch]);
            HilbertMag = HilbertMagPhaseMean[0];
            HilbertPhase = HilbertMagPhaseMean[1];
            delete[] HilbertMagPhaseMean;
            
            HilbertPhase += fieldPhase;
            for (int i=0; i < fReceiverHandler.GetFilterSize(); i++)  // populate filter with field.
            {
                HilbertPhase += 2.*LMCConst::Pi()*fieldfrequency*fReceiverHandler.GetFilterResolution();
                PatchFIRBuffer[channel*fNElementsPerStrip+patch].push_back(HilbertMag*cos(HilbertPhase));
                PatchFIRBuffer[channel*fNElementsPerStrip+patch].pop_front();
            }
            
            double convolution=fReceiverHandler.ConvolveWithFIRFilter(PatchFIRBuffer[channel*fNElementsPerStrip+patch]);
            
            PatchFIRBuffer[channel*fNElementsPerStrip+patch].shrink_to_fit();  // memory deallocation.
            return convolution;
        }
        else return 0.;
        
    }
    
    // Should be moved to somewhere more general
    double TurnstileSignalGenerator::RotateZ(int component, double angle, double x, double y)
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


    bool TurnstileSignalGenerator::InitializePowerCombining()
    {
    	fPowerCombiner.SetSMatrixParameters(fNElementsPerStrip);
    	if (!fPowerCombiner.SetVoltageDampingFactors(fNElementsPerStrip, fElementSpacing) )
    	{
    		return false;
    	}
    	else
    	{
    		return true;
    	}
    	return true;

    }

    
    bool TurnstileSignalGenerator::InitializeElementArray()
    {
        if(!fReceiverHandler.ReadHFSSFile())
        {
            exit(-1);
        }
        const unsigned nChannels = fNChannels;
        const int nReceivers = fNElementsPerStrip;
        const double patchSpacingZ = fElementSpacing;
        const double patchRadius = fArrayRadius;
        double zPosition;
        double theta;
        const double dThetaArray = 2. * LMCConst::Pi() / nChannels; //Divide the circle into nChannels
        const double dRotateVoltages = 0.;  // set to zero to not rotate patch polarities.
        
        Receiver* modelElement = new Receiver;
        modelElement = fPowerCombiner.ChooseElement();

        
        allRxChannels.resize(nChannels);
        
        for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
        {
            theta = channelIndex * dThetaArray;
            
            for(int receiverIndex = 0; receiverIndex < nReceivers; ++receiverIndex)
            {
                zPosition =  (receiverIndex - (nReceivers - 1.) /2.) * patchSpacingZ;
                
                modelElement->SetCenterPosition({patchRadius * cos(theta) , patchRadius * sin(theta) , zPosition });
                modelElement->SetPolarizationDirection({RotateZ(0, dRotateVoltages*channelIndex, sin(theta), -cos(theta)), RotateZ(1, dRotateVoltages*channelIndex, sin(theta), -cos(theta)), 0.});
                
                modelElement->SetNormalDirection({-cos(theta), -sin(theta), 0.}); //Say normals point inwards
                allRxChannels[channelIndex].AddReceiver(modelElement);
            }
        }
        return true;
    }
    
    bool TurnstileSignalGenerator::DoGenerate( Signal* aSignal )
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }
    
    
    bool TurnstileSignalGenerator::DoGenerateTime( Signal* aSignal )
    {

    	if(!InitializeElementArray())
    	{
    		LERROR(lmclog,"Error initilizing Patch Array");
    		exit(-1);
    	}

        if (!InitializePowerCombining() )
        {
        	LERROR(lmclog,"Error configuring Power Combining");
            exit(-1);
        }



        const unsigned nchannels = fNChannels;
        const unsigned nelements = fNElementsPerStrip;
        
        double LO_phase = 0.;
        double field_phase = 0.;
        double VoltageSample = 0.;
        unsigned nfilterbins = fReceiverHandler.GetFilterSize();
        unsigned nfieldbufferbins = fFieldBufferSize;
        double dtfilter = fReceiverHandler.GetFilterResolution();
        
        if(!fAntennaSignalTransmitter.InitializeTransmitter())
        {
            exit(-1);
        }
        
        InitializeBuffers(nfilterbins, nfieldbufferbins);
        
        double initialPhaseDelay=fAntennaSignalTransmitter.GetInitialPhaseDelay();
        for( unsigned index = 0; index < aSignal->TimeSize()*aSignal->DecimationFactor(); ++index )
        {
//            double fieldValue=fAntennaSignalTransmitter.GenerateSignal(aSignal,fAcquisitionRate);
//            double fieldValue=fAntennaSignalTransmitter.GetEFieldAtPoint(currentElement->GetPosition().GetX(), currentElement->GetPosition().GetY(), currentElement->GetPosition().GetZ()).GetZ();
            double antennaPositionX=fAntennaSignalTransmitter.GetAntennaPosition().GetX();
            double antennaPositionY=fAntennaSignalTransmitter.GetAntennaPosition().GetY();
            double antennaPositionZ=fAntennaSignalTransmitter.GetAntennaPosition().GetZ();
            LO_phase += 2.*LMCConst::Pi()*fLO_frequency/aSignal->DecimationFactor()/(fAcquisitionRate*1.e6);
            for (unsigned ch = 0; ch < nchannels; ++ch)
            {
                for (unsigned element = 0; element < nelements; ++element)
                {
                    Receiver *currentElement;
                    currentElement = allRxChannels[ch][element];
                    double* FieldSolution = fAntennaSignalTransmitter.GetEFieldCoPol(currentElement->GetPosition(), currentElement->GetPolarizationDirection());
                    double relativePatchPosX=currentElement->GetPosition().GetX() - antennaPositionX;
                    double relativePatchPosY=currentElement->GetPosition().GetY() - antennaPositionY;
                    double relativePatchPosZ=currentElement->GetPosition().GetZ() - antennaPositionZ;
                    double patchAntennaDistance = sqrt(relativePatchPosX*relativePatchPosX+relativePatchPosY*relativePatchPosY+relativePatchPosZ*relativePatchPosZ);
                    double field_phase=initialPhaseDelay+2.*LMCConst::Pi()*(patchAntennaDistance/LMCConst::C())*fRF_frequency;
                    FillBuffers(aSignal, FieldSolution[0], field_phase, LO_phase, index, ch, element);
                    VoltageSample = GetVoltageFromField(ch, element, field_phase)*GetAOIFactor(currentElement->GetPosition()-fAntennaSignalTransmitter.GetAntennaPosition(),currentElement->GetPosition())/patchAntennaDistance;;
     	            fPowerCombiner.AddOneVoltageToStripSum(aSignal, VoltageSample, LO_phase, element, IndexBuffer[ch*fNElementsPerStrip+element].front());
                    PopBuffers(ch, element);
                }  // patch
            }  // channel
            if ( index%fSwapFrequency == 0 ) CleanupBuffers();  // release memory
        }  // index
        return true;
    }
    
    void TurnstileSignalGenerator::FillBuffers(Signal* aSignal, double FieldValue, double FieldPhase, double LOPhase, unsigned index, unsigned channel, unsigned patch)
    {
        
        EFieldBuffer[channel*fNElementsPerStrip+patch].push_back(FieldValue);
        EPhaseBuffer[channel*fNElementsPerStrip+patch].push_back(FieldPhase);
        EAmplitudeBuffer[channel*fNElementsPerStrip+patch].push_back(fAmplitude);
        EFrequencyBuffer[channel*fNElementsPerStrip+patch].push_back(fRF_frequency);
        LOPhaseBuffer[channel*fNElementsPerStrip+patch].push_back(LOPhase);
        IndexBuffer[channel*fNElementsPerStrip+patch].push_back(channel*aSignal->TimeSize()*aSignal->DecimationFactor() + index);
        
    }
    
    
    void TurnstileSignalGenerator::PopBuffers(unsigned channel, unsigned patch)
    {
        EFieldBuffer[channel*fNElementsPerStrip+patch].pop_front();
        EPhaseBuffer[channel*fNElementsPerStrip+patch].pop_front();
        EAmplitudeBuffer[channel*fNElementsPerStrip+patch].pop_front();
        EFrequencyBuffer[channel*fNElementsPerStrip+patch].pop_front();
        LOPhaseBuffer[channel*fNElementsPerStrip+patch].pop_front();
        IndexBuffer[channel*fNElementsPerStrip+patch].pop_front();
        
        EFieldBuffer[channel*fNElementsPerStrip+patch].shrink_to_fit();
        EPhaseBuffer[channel*fNElementsPerStrip+patch].shrink_to_fit();
        EAmplitudeBuffer[channel*fNElementsPerStrip+patch].shrink_to_fit();
        EFrequencyBuffer[channel*fNElementsPerStrip+patch].shrink_to_fit();
        LOPhaseBuffer[channel*fNElementsPerStrip+patch].shrink_to_fit();
        IndexBuffer[channel*fNElementsPerStrip+patch].shrink_to_fit();
        // PTS: Seg faults when shrink-to-fit used on ConvolutionTimeBuffer, removed for now. Need to revisit what the problem is
        //ConvolutionTimeBuffer[channel+fNElementsPerStrip+patch].shrink_to_fit();
    }
    
    
    
    
    void TurnstileSignalGenerator::CleanupBuffers()
    {
        FieldBuffer aFieldBuffer;
        EFieldBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
        EPhaseBuffer = aFieldBuffer.CleanupBuffer(EPhaseBuffer);
        EAmplitudeBuffer = aFieldBuffer.CleanupBuffer(EAmplitudeBuffer);
        EFrequencyBuffer = aFieldBuffer.CleanupBuffer(EFrequencyBuffer);
        LOPhaseBuffer = aFieldBuffer.CleanupBuffer(LOPhaseBuffer);
        IndexBuffer = aFieldBuffer.CleanupBuffer(IndexBuffer);
        
    }
    
    void TurnstileSignalGenerator::InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize)
    {
        
        FieldBuffer aFieldBuffer;
        
        EFieldBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, fieldbuffersize);
        EPhaseBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, fieldbuffersize);
        EAmplitudeBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, fieldbuffersize);
        EFrequencyBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, fieldbuffersize);
        LOPhaseBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, fieldbuffersize);
        IndexBuffer = aFieldBuffer.InitializeUnsignedBuffer(fNChannels, fNElementsPerStrip, fieldbuffersize);
        PatchFIRBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNElementsPerStrip, filterbuffersize);
    }
    
    
    bool TurnstileSignalGenerator::DoGenerateFreq( Signal* aSignal )
    {
        return true;
    }
    
} /* namespace locust */

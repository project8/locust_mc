
/*
 * LMCPlaneWaveSignalGenerator.cc
 *
 *  Created on: Feb 19, 2018
 *      Author: a. b. telles
 */

#include "LMCPlaneWaveSignalGenerator.hh"

#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <math.h>       /* sin */

#include "LMCGlobalsDeclaration.hh"
#include "LMCDigitizer.hh"


namespace locust
{
	LOGGER( lmclog, "PlaneWaveSignalGenerator" );
	MT_REGISTER_GENERATOR(PlaneWaveSignalGenerator, "planewave-signal");
	PlaneWaveSignalGenerator::PlaneWaveSignalGenerator( const std::string& aName ) :
    	Generator( aName ),
		fLO_Frequency( 0.),
		fphiLO(0.),
		fRF_Frequency( 0.),
		fArrayRadius( 0. ),
		fNPatchesPerStrip( 0. ),
		fPatchSpacing( 0. ),
		fFieldBufferSize( 50 ),
		fAOI( 0.),
		fAmplitude( 0.),
		fSwapFrequency( 1000 )

	{
		fRequiredSignalState = Signal::kTime;
	}

	PlaneWaveSignalGenerator::~PlaneWaveSignalGenerator()
	{
	}

	bool PlaneWaveSignalGenerator::Configure( const scarab::param_node& aParam )
	{
      
		if(!fTFReceiverHandler.Configure(aParam))
		{
			LERROR(lmclog,"Error configuring receiver TFHandler class");
		}

		if(!fPowerCombiner.Configure(aParam))
		{
			LERROR(lmclog,"Error configuring PowerCombiner class");
		}

		if( aParam.has( "buffer-size" ) )
		{
			fFieldBufferSize = aParam["buffer-size"]().as_int();
			fHilbertTransform.SetBufferSize(aParam["buffer-size"]().as_int());
		}

		if(!fHilbertTransform.Configure(aParam))
		{
			LERROR(lmclog,"Error configuring HilbertTransform class");
		}

		if( aParam.has( "planewave-frequency" ) )
		{
			SetPlaneWaveFrequency( aParam.get_value< double >( "planewave-frequency", fRF_Frequency ));
		}

		if( aParam.has( "lo-frequency" ) )
		{
			SetLOFrequency( aParam.get_value< double >( "lo-frequency", fLO_Frequency ));
		}

		if( aParam.has( "array-radius" ) )
		{
			SetArrayRadius( aParam.get_value< double >( "array-radius", fArrayRadius ));
		}

		if( aParam.has( "npatches-per-strip" ) )
		{
			SetNPatchesPerStrip( aParam.get_value< int >( "npatches-per-strip", fNPatchesPerStrip ));
		}
		if( aParam.has( "patch-spacing" ) )
		{
			SetPatchSpacing( aParam.get_value< double >( "patch-spacing", fPatchSpacing ) );
		}
		if( aParam.has( "AOI" ) )
		{
			SetAOI( aParam.get_value< double >( "AOI", fAOI ));
		}
		if( aParam.has( "amplitude" ) )
		{
			SetAmplitude( aParam.get_value< double >( "amplitude", fAmplitude) );
		}
		return true;
	}

    double PlaneWaveSignalGenerator::GetPlaneWaveFrequency() const
    {
    	return fRF_Frequency;
    }

    void PlaneWaveSignalGenerator::SetPlaneWaveFrequency( double aPlaneWaveFrequency )
    {
        fRF_Frequency = aPlaneWaveFrequency;
        return;
    }

    double PlaneWaveSignalGenerator::GetLOFrequency() const
    {
    	return fLO_Frequency;
    }

    void PlaneWaveSignalGenerator::SetLOFrequency( double aLOFrequency )
    {
    	fLO_Frequency = aLOFrequency;
    	return;
    }

    double PlaneWaveSignalGenerator::GetArrayRadius() const
    {
    	return fArrayRadius;
    }

    void PlaneWaveSignalGenerator::SetArrayRadius( double aArrayRadius )
    {
    	fArrayRadius = aArrayRadius;
        return;
    }

    int PlaneWaveSignalGenerator::GetNPatchesPerStrip() const
    {
    	return fNPatchesPerStrip;
    }

    void PlaneWaveSignalGenerator::SetNPatchesPerStrip( int aNPatchesPerStrip )
    {
    	fNPatchesPerStrip = aNPatchesPerStrip;
    	return;
    }

    double PlaneWaveSignalGenerator::GetPatchSpacing() const
    {
    	return fPatchSpacing;
    }

    void PlaneWaveSignalGenerator::SetPatchSpacing( double aPatchSpacing )
    {
    	fPatchSpacing = aPatchSpacing;
        return;
    }

    double PlaneWaveSignalGenerator::GetAOI() const
    {
    	return fAOI;
    }

    void PlaneWaveSignalGenerator::SetAOI( double aAOI )
    {
    	fAOI = aAOI*(2*LMCConst::Pi()/360); //convert to radians
        return;
    }
    double PlaneWaveSignalGenerator::GetAmplitude() const
    {
    	return fAmplitude;
    }

    void PlaneWaveSignalGenerator::SetAmplitude( double aAmplitude )
    {
    	fAmplitude = aAmplitude;
    	return;
    }


    void PlaneWaveSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
    	aVisitor->Visit( this );
    	return;
    }

    double PlaneWaveSignalGenerator::GetAOIFactor(double AOI, LMCThreeVector PatchNormalVector)
    {
    	LMCThreeVector IncidentKVector;
    	IncidentKVector.SetComponents(cos(AOI), 0.0, sin(AOI));
    	double AOIFactor = fabs(IncidentKVector.Dot(PatchNormalVector));
    	return AOIFactor;
    }


    double PlaneWaveSignalGenerator::GetPWPhaseDelayAtPatch(int z_index)
    {
    	double phasedelay = 0.;
    	if(fAOI >= 0)
    	{
    		phasedelay = 2*LMCConst::Pi()*z_index*fPatchSpacing*sin(fAOI)*fRF_Frequency/LMCConst::C();
    	}
    	else
    	{
    		phasedelay = (fNPatchesPerStrip - z_index)*2*LMCConst::Pi()*fPatchSpacing*sin(fAOI)*fRF_Frequency/LMCConst::C();
    	}
    	return phasedelay;
    }
  
    double PlaneWaveSignalGenerator::GetPatchFIRSample(double dottedamp, double startphase, int patchIndex)
    {
   
//    	double* generatedpoints = new double [nfilterbins];
    	std::deque<double> generatedpoints;
    	int nfilterbins = fTFReceiverHandler.GetFilterSize();
    	double dtfilter = fTFReceiverHandler.GetFilterResolution();
    	//double phase = startphase;
    	double phase = startphase + GetPWPhaseDelayAtPatch(patchIndex);
    	double amp = dottedamp;
    
    	for(int i=0; i < nfilterbins; i++)
    	{
    		generatedpoints.push_back(amp*cos(phase));
    		phase += 2*LMCConst::Pi()*dtfilter*fRF_Frequency;

    		// TEST PRINT STATEMENT
	//		printf("genpoints %d is %g, amp is %g\n", i, generatedpoints[i], amp); getchar();
    	}

    	double convolution=fTFReceiverHandler.ConvolveWithFIRFilter(generatedpoints);
      
    	generatedpoints.shrink_to_fit();  // memory deallocation.

    	return convolution;
    }
  
    double* PlaneWaveSignalGenerator::GetHilbertMagPhase(unsigned bufferIndex)
    {
    	double* magphase = new double[2];
    	magphase[0] = 0.;
    	magphase[1] = 0.;
    
    	if (fabs(PWValueBuffer[bufferIndex].front()) > 0.)
    	{
    		HilbertTransform aHilbertTransform;
    		double* HilbertMagPhaseMean = new double[3];
	
    		HilbertMagPhaseMean = aHilbertTransform.GetMagPhaseMean(PWValueBuffer[bufferIndex], PWFreqBuffer[bufferIndex]);
    		magphase[0] = HilbertMagPhaseMean[0];
    		magphase[1] = HilbertMagPhaseMean[1];
    		delete[] HilbertMagPhaseMean;
    	}
    	return magphase;
    }
  

    void PlaneWaveSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal)
    {
    	unsigned bufferIndex = 0;
    	const int signalSize = aSignal->TimeSize();
    	const double timeSampleSize = 1./(1.e6 * fAcquisitionRate * aSignal->DecimationFactor());
    	unsigned sampleIndex = 0;
    	double fieldamp = 0.;
    	double fieldphase = 0.;
    	double fieldvalue = 0.;
    	double* hilbertmagphase = new double[2];
    	fphiLO += 2. * LMCConst::Pi() * fLO_Frequency * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

    	for(int channelIndex = 0; channelIndex < allChannels.size(); ++channelIndex)
    	{
    		for(int patchIndex = 0; patchIndex < allChannels[channelIndex].size(); ++patchIndex)
    		{
    			sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;
    			bufferIndex = channelIndex*fNPatchesPerStrip+patchIndex;
	    
    			PatchAntenna *currentPatch;
    			currentPatch = &allChannels[channelIndex][patchIndex];

    			fieldamp = fAmplitude*GetAOIFactor(fAOI, currentPatch->GetNormalDirection());
    			fieldphase = PWPhaseBuffer[bufferIndex].back();
    			fieldphase += 2. * LMCConst::Pi() * fRF_Frequency * timeSampleSize;
	    // 		fieldphase += GetPWPhaseDelayAtPatch(patchIndex); this really does not work here for some reason. need to investigate further.
    			fieldvalue = fieldamp*cos(fieldphase);

    			FillBuffers(bufferIndex, sampleIndex, fieldphase, fieldvalue);
    			PopBuffers(bufferIndex);

    			hilbertmagphase = GetHilbertMagPhase(bufferIndex);
    			int hilbertbuffermargin = fHilbertTransform.GetBufferMargin();

    			PatchVoltageBuffer[bufferIndex].emplace(PatchVoltageBuffer[bufferIndex].begin()+hilbertbuffermargin+1, GetPatchFIRSample(hilbertmagphase[0], hilbertmagphase[1], patchIndex));
    			PatchVoltageBuffer[bufferIndex].pop_front();
    			PatchVoltageBuffer[bufferIndex].shrink_to_fit();
    			fPowerCombiner.AddOneVoltageToStripSum(aSignal, PatchVoltageBuffer[bufferIndex].front(), fphiLO, patchIndex, SampleIndexBuffer[bufferIndex].front());
	   
    			// TEST PRINT STATEMENTS
/*
	      		printf("Channel is %d\n", channelIndex);
	      	  	printf("Patch is %d\n", patchIndex);
	      	  	printf("Digitizer Sample is %d\n", index);

	      	    printf("fieldamp is %f\n", fieldamp);
	      	  	printf("fieldphase is %f\n", fieldphase);
	      	  	printf("fieldvalue is %f\n", fieldvalue);
	      	  	printf("hilbertmagphase[0] is %g\n", hilbertmagphase[0]);
	      	  	printf("hilbertmagphase[1] is %g\n", hilbertmagphase[1]);
	    
	      	  	printf("SampleIndexBuffer[%d] is %u\n", bufferIndex, SampleIndexBuffer[bufferIndex].front());
	      	  	printf("LOPhaseBuffer[%d] is %f\n", bufferIndex, LOPhaseBuffer[bufferIndex].front());
	      	  	printf("PWFreqBuffer[%d] is %f\n", bufferIndex, PWFreqBuffer[bufferIndex].front());
	      	  	printf("PWPhaseBuffer[%d] is %f\n", bufferIndex, PWPhaseBuffer[bufferIndex].front());
	      	  	printf("PWValueBuffer[%d] is %f\n", bufferIndex, PWValueBuffer[bufferIndex].front());
	      	  	printf("PatchVoltageBuffer[%d] is %e\n", bufferIndex, PatchVoltageBuffer[bufferIndex].front());
	      	  	printf("Resulting VI[%d] is %e\n", sampleIndex, aSignal->LongSignalTimeComplex()[sampleIndex][0]);
	    
	      	  	getchar();
*/

    			//text file for hilbert transform testing.
	    /*
	      	  	std::ofstream hilbertfile;
	      	  	hilbertfile.open("hilbertfile.txt", std::fstream::app);
	      	  	if(patchIndex == 0){
	      		hilbertfile << PWValueBuffer[bufferIndex].front();
	      	  	hilbertfile << ", ";
	      	  	hilbertfile << PWPhaseBuffer[bufferIndex].front();
	      	  	hilbertfile << ", ";
	      	  	hilbertfile << hilbertmagphase[0];
	      	  	hilbertfile << ", ";
	      	  	hilbertfile << hilbertmagphase[1];
	      	  	hilbertfile << "\n";
	      	  	hilbertfile.close();

	    */
	    
    		}  // patch
    	} // channel
        if ( index%fSwapFrequency == 0 ) CleanupBuffers();  // release memory

    }
    
  

    void PlaneWaveSignalGenerator::FillBuffers(unsigned bufferIndex, int digitizerIndex, double pwphase, double pwval)
    {
    	SampleIndexBuffer[bufferIndex].push_back(digitizerIndex);
    	PWFreqBuffer[bufferIndex].push_back(fRF_Frequency);
    	PWPhaseBuffer[bufferIndex].push_back(pwphase);
    	PWValueBuffer[bufferIndex].push_back(pwval);
    }

    void PlaneWaveSignalGenerator::PopBuffers(unsigned bufferIndex)
    {
    	SampleIndexBuffer[bufferIndex].pop_front();
    	PWFreqBuffer[bufferIndex].pop_front();
    	PWPhaseBuffer[bufferIndex].pop_front();
    	PWValueBuffer[bufferIndex].pop_front();
    }

    void PlaneWaveSignalGenerator::CleanupBuffers()
    {
    	FieldBuffer aFieldBuffer;
    	PWValueBuffer = aFieldBuffer.CleanupBuffer(PWValueBuffer);
    	PWFreqBuffer = aFieldBuffer.CleanupBuffer(PWFreqBuffer);
    	PWPhaseBuffer = aFieldBuffer.CleanupBuffer(PWPhaseBuffer);
    	SampleIndexBuffer = aFieldBuffer.CleanupBuffer(SampleIndexBuffer);
    }

  
    void PlaneWaveSignalGenerator::InitializeBuffers()
    {
    	const unsigned nchannels = fNChannels;
    	const int nReceivers = fNPatchesPerStrip;
    
    	FieldBuffer aFieldBuffer;

    	SampleIndexBuffer = aFieldBuffer.InitializeUnsignedBuffer(nchannels, nReceivers, fFieldBufferSize);
    	PWFreqBuffer = aFieldBuffer.InitializeBuffer(nchannels, nReceivers, fFieldBufferSize);
    	PWValueBuffer = aFieldBuffer.InitializeBuffer(nchannels, nReceivers, fFieldBufferSize);
    	PWPhaseBuffer = aFieldBuffer.InitializeBuffer(nchannels, nReceivers, fFieldBufferSize);
    	PatchVoltageBuffer = aFieldBuffer.InitializeBuffer(nchannels, nReceivers, fFieldBufferSize);
    
    }



    bool PlaneWaveSignalGenerator::InitializePowerCombining()
    {
    	fPowerCombiner.SetSMatrixParameters(fNPatchesPerStrip);
    	fPowerCombiner.SetVoltageDampingFactors(fNPatchesPerStrip);

    	return true;

    }



    bool PlaneWaveSignalGenerator::InitializePatchArray()
    {
    	if(!fTFReceiverHandler.ReadHFSSFile())
    	{
    		return false;
    	}
    	const unsigned nChannels = fNChannels;
    	const int nReceivers = fNPatchesPerStrip;

    	const double patchSpacingZ = fPatchSpacing;
    	const double patchRadius = fArrayRadius;
    	double zPosition;
    	double theta;
    	const double dThetaArray = 2. * LMCConst::Pi() / nChannels; //Divide the circle into nChannels

    	PatchAntenna modelPatch;

    	allChannels.resize(nChannels);

    	for(int channelIndex = 0; channelIndex < nChannels; ++channelIndex)
    	{
    		theta = channelIndex * dThetaArray;

    		for(int receiverIndex = 0; receiverIndex < nReceivers; ++receiverIndex)
    		{
    			zPosition =  (receiverIndex - (nReceivers - 1.) /2.) * patchSpacingZ;

                if (fPowerCombiner.GetPowerCombiner() == 7)  // single patch
                {
                	zPosition = 0.;
                }

                modelPatch.SetCenterPosition({patchRadius * cos(theta) , patchRadius * sin(theta) , zPosition });
    			modelPatch.SetPolarizationDirection({sin(theta), -cos(theta), 0.});
    			modelPatch.SetNormalDirection({-cos(theta), -sin(theta), 0.}); //Say normals point inwards
    			allChannels[channelIndex].AddReceiver(modelPatch);
    		}
    	}
    	return true;
    }



    bool PlaneWaveSignalGenerator::DoGenerate( Signal* aSignal )
    {

    	if(!InitializePatchArray())
	{
	    LERROR(lmclog,"Error initilizing Patch Array");
	    exit(-1);
	}
    	InitializePowerCombining();

    	int nfilterbins = fTFReceiverHandler.GetFilterSize();
 
    	InitializeBuffers();

    	// text file for VI for testing.
    	// note: this is pre-low pass filter. better to put this into LMCDecimateSignalGenerator.cc
    	//   std::ofstream voltagefile;
    	//    voltagefile.open("voltagefile.txt");

    	//n samples for event spacing.
    	int PreEventCounter = 0;
    	const int NPreEventSamples = 150000;
    	PreEventCounter = NPreEventSamples; // jump past wait time.

    	for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
    	{
    		DriveAntenna(PreEventCounter, index, aSignal);
    		/*
	  	  	voltagefile << index;
	  	  	voltagefile << "\n";
	  	  	voltagefile << aSignal->LongSignalTimeComplex()[index][0];
	  	  	voltagefile << "\n";
	  	  	voltagefile << aSignal->LongSignalTimeComplex()[index][1];
	  	  	voltagefile << "\n";
    		 */
    	}
    	//  voltagefile.close();


    	return true;
    }

} /* namespace locust */

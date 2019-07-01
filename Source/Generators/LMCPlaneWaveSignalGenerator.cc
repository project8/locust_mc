/*
 * LMCPlaneWaveSignalGenerator.cc
 *
 *  Created on: Feb 19, 2018
 *      Author: pslocum
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
    fRF_Frequency( 0.),
    fArrayRadius( 0. ),
    fPhaseDelay( 0 ),
    fVoltageDamping( 0 ),
    fNPatchesPerStrip( 0. ),
    fPatchSpacing( 0. ),
    fPowerCombiner( 0 ),
    fAOI( 0.),
    fPatchFIRfilter( 0.),
    gpatchfilter_filename("blank.txt"),
    fPatchFIRfilter_resolution( 0. ),
    fAmplitude( 0.)
  {
    fRequiredSignalState = Signal::kTime;
  }

  PlaneWaveSignalGenerator::~PlaneWaveSignalGenerator()
  {
  }

  bool PlaneWaveSignalGenerator::Configure( const scarab::param_node* aParam )
  {
    if( aParam == NULL) return true;

    if( aParam->has( "phase-delay" ) )
      {
	fPhaseDelay = aParam->get_value< bool >( "phase-delay" );
      }

    if( aParam->has( "voltage-damping" ) )
      {
	fVoltageDamping = aParam->get_value< bool >( "voltage-damping" );
      }

    if( aParam->has( "planewave-frequency" ) )
      {
	fRF_Frequency = aParam->get_value< double >( "planewave-frequency" );
      }

    if( aParam->has( "lo-frequency" ) )
      {
	fLO_Frequency = aParam->get_value< double >( "lo-frequency" );
      }
    if( aParam->has( "array-radius" ) )
      {
	fArrayRadius = aParam->get_value< double >( "array-radius" );
      }
    if( aParam->has( "npatches-per-strip" ) )
      {
	fNPatchesPerStrip = aParam->get_value< int >( "npatches-per-strip" );
      }
    if( aParam->has( "patch-spacing" ) )
      {
	fPatchSpacing = aParam->get_value< double >( "patch-spacing" );
      }
    if( aParam->has( "feed" ) )
      {
	if (aParam->get_value< std::string >( "feed" ) == "corporate")
	  fPowerCombiner = 0;  // default
	else if (aParam->get_value< std::string >( "feed" ) == "series")
	  fPowerCombiner = 1;
	else if (aParam->get_value< std::string >( "feed" ) == "one-quarter")
	  fPowerCombiner = 2;
	else if (aParam->get_value< std::string >( "feed" ) == "seven-eighths")
	  fPowerCombiner = 3;
	else if (aParam->get_value< std::string >( "feed") == "nine-sixteenths")
	  fPowerCombiner = 4;
	else
	  fPowerCombiner = 0;  // default
      }
    if( aParam->has( "AOI" ) )
      {
	fAOI = aParam->get_value< double >( "AOI" );
	fAOI *= (2*LMCConst::Pi()/360); //convert to radians
      }
    if ( aParam->has(  "patch-filter"  ) )
      {
	fPatchFIRfilter = aParam->get_value< bool >( "patch-filter" );
      }
    if( aParam->has( "patch-filter-filename" ) )
      {
	gpatchfilter_filename = aParam->get_value< std::string >( "patch-filter-filename" );
      }
    if( aParam->has( "patch-filter-resolution" ) )
      {
	fPatchFIRfilter_resolution = aParam->get_value< double >( "patch-filter-resolution" );
      }
    if( aParam->has( "amplitude" ) )
      {
	fAmplitude = aParam->get_value< double >( "amplitude" );
      }

    return true;
  }

  void PlaneWaveSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
  {
    aVisitor->Visit( this );
    return;
  }

  double PlaneWaveSignalGenerator::GetAOIFactor(double AOI, LMCThreeVector PatchNormalVector){
    LMCThreeVector IncidentKVector;
    IncidentKVector.SetComponents(cos(AOI), 0.0, sin(AOI));
    double AOIFactor = fabs(IncidentKVector.Dot(PatchNormalVector));
    return AOIFactor;
  }


  double PlaneWaveSignalGenerator::GetVoltageAmpFromPlaneWave(int z_index)
  {
    double AntennaFactor = 0.;
    double amplitude = 0.;
    if(z_index == 0 || z_index == fNPatchesPerStrip)
      {
	AntennaFactor = 1./430.;
      }
    else
      {
	AntennaFactor = 1./460.;
      }

    // S = epsilon0 c E0^2 / 2.  // power/area
    //  0.6e-21 W/Hz * 24.e3 Hz / (0.00375*0.002916) = S = 1.3e-12 W/m^2
    // We should detect 0.6e-21 W/Hz * 24.e3 Hz in Katydid.
    // E0 = sqrt(2.*S/epsilon0/c)
    // effective patch area 0.00004583662 m^2

    double S = 0.6e-21*24.e3/(0.00004271);  // W/m^2, effective aperture.
    double E0 = sqrt(2.*S/(LMCConst::EpsNull() * LMCConst::C()));
    //    double E0 = 1.0; // V/m, test case
    amplitude = E0*AntennaFactor;  // volts
      
    // printf("amplitude is %g\n", amplitude); getchar();
    return amplitude;
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
  
  void PlaneWaveSignalGenerator::ProcessFIRFilter(int nskips)
  {

    FILE *fp;
    double filter;
    double index;
    fp = fopen(gpatchfilter_filename.c_str(),"r");
    int count = 0;


    for (int i=0; i<sizeof(FIR_array)/sizeof(FIR_array[0]); i++)
      FIR_array[i] = -99.;


    while (!feof(fp))
      {
        fscanf(fp, "%lf %lf\n", &index, &filter);
        if (count%nskips==0) FIR_array[count/nskips] = filter;
        count += 1;
      }

    fclose(fp);

  }

  int PlaneWaveSignalGenerator::GetNFilterBins()
  {
    int nbins = 0;
    for (int i=0; i<sizeof(FIR_array)/sizeof(FIR_array[0]); i++)
      {
	if (FIR_array[i]>0.) nbins += 1;

	// TEST
	// printf("FIR_array[%d] is %g\n", i, FIR_array[i]); getchar();
	
      }    
    return nbins;
  }
  

  double PlaneWaveSignalGenerator::GetPatchFIRSample(double dottedamp, double startphase, int patchIndex)
  {   
   
    double* generatedpoints = new double [nfilterbins];
    double dtfilter = fPatchFIRfilter_resolution;
    double phase = startphase + GetPWPhaseDelayAtPatch(patchIndex);
    double amp = dottedamp;
    
    for(int i=0; i < nfilterbins; i++)
      {
	generatedpoints[i] = amp*cos(phase);
	phase += 2*LMCConst::Pi()*dtfilter*fRF_Frequency;

	// TEST PRINT STATEMENT
	//	printf("genpoints %d is %g", i, generatedpoints[i]); getchar();
      }

    double total = 0.;
    for(int j=0; j < nfilterbins; j++)
      {
	total += generatedpoints[j]*FIR_array[j];
      }
    delete[] generatedpoints;
    return total;
      
  }

  
  double* PlaneWaveSignalGenerator::GetHilbertMagPhase(unsigned bufferIndex)
  {
    double* magphase = new double[2];
    magphase[0] = 0.;
    magphase[1] = 0.;
    
    if (fabs(PWMagBuffer[bufferIndex].front()) > 0.)
      {
	double fFieldBufferMargin = 50;
	HilbertTransform aHilbertTransform;
	double* HilbertMagPhaseMean = new double[3];
	
	HilbertMagPhaseMean = aHilbertTransform.GetMagPhaseMean(PWMagBuffer[bufferIndex], PWFreqBuffer[bufferIndex], fFieldBufferMargin, 1.e6*fAcquisitionRate*10);
	magphase[0] = HilbertMagPhaseMean[0];
	magphase[1] = HilbertMagPhaseMean[1];
	delete[] HilbertMagPhaseMean;
    
      }
    return magphase;
  }
  

  // z-index ranges from 0 to npatches-per-strip-1.
  void PlaneWaveSignalGenerator::AddOnePatchVoltageToStripSum(Signal* aSignal, unsigned bufferIndex, int patchIndex)
  {
    unsigned sampleIndex = SampleIndexBuffer[bufferIndex].front();
    double phi_LO = LOPhaseBuffer[bufferIndex].front();
    double VoltagePhase = 0.; // need to fix this later
    double VoltageAmplitude = PatchVoltageBuffer[bufferIndex].front();
    unsigned z_index = patchIndex;
    double DopplerFrequency = fRF_Frequency;
    
    PowerCombiner aPowerCombiner;
    
    if (fPowerCombiner == 0 ) //corporate feed, for testing
      {
	if (fVoltageDamping)
	  {
	    VoltageAmplitude *= aPowerCombiner.GetCorporateVoltageDamping();
	  }
      }

    if (fPowerCombiner == 1)  // series feed
      {
	if (fPhaseDelay) // parameter from json
	  {
	    VoltagePhase += aPowerCombiner.GetSeriesPhaseDelay(z_index, DopplerFrequency, fPatchSpacing);
	  }
	if (fVoltageDamping)  // parameter from json
	  {
	    VoltageAmplitude *= aPowerCombiner.GetSeriesVoltageDamping(z_index);
	  }
      }

    if (fPowerCombiner == 2) // one-quarter power combining, center fed strip
      {
	if (fPhaseDelay)
	  {
	    VoltagePhase += aPowerCombiner.GetCenterFedPhaseDelay(fNPatchesPerStrip, z_index, DopplerFrequency, fPatchSpacing);
	  }
	if (fVoltageDamping)
	  {
	    VoltageAmplitude *= aPowerCombiner.GetOneQuarterVoltageDamping(fNPatchesPerStrip, z_index);
	  }
      }

    if (fPowerCombiner == 3) // seven-eighths power combining, center fed strip
      {
	if (fPhaseDelay)
	  {
	    VoltagePhase += aPowerCombiner.GetCenterFedPhaseDelay(fNPatchesPerStrip, z_index, DopplerFrequency, fPatchSpacing);
	  }
	if (fVoltageDamping)
	  {
	    VoltageAmplitude *= aPowerCombiner.GetSevenEighthsVoltageDamping(fNPatchesPerStrip, z_index);
	  }
      }

    if (fPowerCombiner == 4) // nine-sixteenths power combining, center fed strip
      {
	if (fPhaseDelay)
	  {
	    VoltagePhase += aPowerCombiner.GetCenterFedPhaseDelay(fNPatchesPerStrip, z_index, DopplerFrequency, fPatchSpacing);
	  }
	if (fVoltageDamping)
	  {
	    VoltageAmplitude *= aPowerCombiner.GetNineSixteenthsVoltageDamping(fNPatchesPerStrip, z_index);
	  }
      }

    // factor of 2 is needed for cosA*cosB = 1/2*(cos(A+B)+cos(A-B)); usually we leave out the 1/2 for e.g. sinusoidal RF.
    
    aSignal->LongSignalTimeComplex()[sampleIndex][0] += VoltageAmplitude * 2. * cos(phi_LO);
    aSignal->LongSignalTimeComplex()[sampleIndex][1] += VoltageAmplitude * 2. * cos(LMCConst::Pi()/2 + phi_LO);

    /*
     
      double VI = 0.;
      double VQ = 0.;
      double phi_at_patch = 0.;

      VI = VoltageAmplitude * 2. * cos(phi_LO);
      VQ = VoltageAmplitude * 2. * cos(LMCConst::Pi()/2 + phi_LO);

     
      phi_at_patch = atan(abs(VQ)/abs(VI));

      if(VQ > 0. && VI < 0.) // second quadrant
      {
      phi_at_patch = LMCConst::Pi()-phi_at_patch;
      }
      else if(VQ < 0. && VI < 0.) // third quadrant
      {
      phi_at_patch = LMCConst::Pi()+phi_at_patch;
      }
      else if(VQ < 0. && VI > 0.) // fourth quadrant
      {
      phi_at_patch = 2*LMCConst::Pi()-phi_at_patch;
      }
     
      if(fPhaseDelay)
      {
	 
      }
    */
    // TEST
    /*
      printf("Voltage Amplitude for patch at %d is %g\n", sampleIndex, VoltageAmplitude);
      printf("summedvoltageamplitude at %d is %g\n", sampleIndex, aSignal->LongSignalTimeComplex()[sampleIndex][0]);
      getchar();
    */

  }



  void* PlaneWaveSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal){
       
    unsigned bufferIndex = 0;
    const int signalSize = aSignal->TimeSize();
    const double timeSampleSize = 1./(1.e6 * fAcquisitionRate * aSignal->DecimationFactor());
    unsigned sampleIndex = 0;
    double fieldamp = 0.;
    double fieldphase = 0.;
    double fieldvalue = 0.;
    double* hilbertmagphase = new double[2];
    double phiLO = 0.;

    for(int channelIndex = 0; channelIndex < allChannels.size(); ++channelIndex)
      {
	for(int patchIndex = 0; patchIndex < allChannels[channelIndex].size(); ++patchIndex)
	  {
	    sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;
	    bufferIndex = channelIndex*fNPatchesPerStrip+patchIndex;

	    phiLO = LOPhaseBuffer[bufferIndex].back();
	    phiLO +=  2. * LMCConst::Pi() * fLO_Frequency * timeSampleSize;
	    
	    PatchAntenna *currentPatch;
	    currentPatch = &allChannels[channelIndex][patchIndex];

	    fieldamp = fAmplitude*GetAOIFactor(fAOI, currentPatch->GetNormalDirection());
	    fieldphase = PWPhaseBuffer[bufferIndex].back();
	    fieldphase += 2. * LMCConst::Pi() * fRF_Frequency * timeSampleSize;
	    fieldvalue = fieldamp*cos(fieldphase);

	    FillBuffers(bufferIndex, sampleIndex, phiLO, fieldphase, fieldvalue);

	    hilbertmagphase = GetHilbertMagPhase(bufferIndex);
	    
	    PatchVoltageBuffer[bufferIndex].push_back(GetPatchFIRSample(hilbertmagphase[0], hilbertmagphase[1], patchIndex));
	    
	    AddOnePatchVoltageToStripSum(aSignal, bufferIndex, patchIndex);
	    	    
	    

	    // TEST PRINT STATEMENTS
	    /*
	      printf("Channel is %d\n", channelIndex);
	      printf("Patch is %d\n", patchIndex);
	      printf("SampleIndexBuffer[%d] is %u\n", bufferIndex, SampleIndexBuffer[bufferIndex].front());
	      printf("LOPhaseBuffer[%d] is %f\n", bufferIndex, LOPhaseBuffer[bufferIndex].front());
	      printf("PWFreqBuffer[%d] is %f\n", bufferIndex, PWFreqBuffer[bufferIndex].front());
	      printf("PWPhaseBuffer[%d] is %f\n", bufferIndex, PWPhaseBuffer[bufferIndex].front());
	      printf("PWMagBuffer[%d] is %f\n", bufferIndex, PWMagBuffer[bufferIndex].front());
	      printf("PatchVoltageBuffer[%d] is %f\n", bufferIndex, PatchVoltageBuffer[bufferIndex].front());
	      printf("Resulting VI[%d] is %f\n", sampleIndex, aSignal->LongSignalTimeComplex()[sampleIndex][0]);
	      //getchar();
	    

	      // TEST HILBERT PRINT STATEMENTS

	      printf("fieldamp is %f\n", fieldamp);
	      printf("fieldphase is %f\n", fieldphase);
	      printf("fieldvalue is %f\n", fieldvalue);
	      printf("hilbertmagphase[0] is %f\n", hilbertmagphase[0]);
	      printf("hilbertmagphase[1] is %f\n", hilbertmagphase[1]);
	      getchar();
	    */

	    //text file for hilbert transform testing.
	    /*
	      std::ofstream hilbertfile;
	      hilbertfile.open("hilbertfile.txt", std::fstream::app);
	      if(patchIndex == 0){
	      hilbertfile << PWMagBuffer[bufferIndex].front();
	      hilbertfile << ", ";
	      hilbertfile << PWPhaseBuffer[bufferIndex].front();
	      hilbertfile << ", ";
	      hilbertfile << hilbertmagphase[0];
	      hilbertfile << ", ";
	      hilbertfile << hilbertmagphase[1];
	      hilbertfile << "\n";
	      hilbertfile.close();
	      }
	    */
	    
	    PopBuffers(bufferIndex);
	  }
      }
  }
    
  

  void PlaneWaveSignalGenerator::FillBuffers(unsigned bufferIndex, int digitizerIndex, double phiLO, double pwphase, double pwmag){

    SampleIndexBuffer[bufferIndex].push_back(digitizerIndex);
    LOPhaseBuffer[bufferIndex].push_back(phiLO);
    PWFreqBuffer[bufferIndex].push_back(fRF_Frequency);
    PWPhaseBuffer[bufferIndex].push_back(pwphase);
    PWMagBuffer[bufferIndex].push_back(pwmag);
  }

  void PlaneWaveSignalGenerator::PopBuffers(unsigned bufferIndex){
    
    SampleIndexBuffer[bufferIndex].pop_front();
    LOPhaseBuffer[bufferIndex].pop_front();
    PWFreqBuffer[bufferIndex].pop_front();
    PWPhaseBuffer[bufferIndex].pop_front();
    PWMagBuffer[bufferIndex].pop_front();
    PatchVoltageBuffer[bufferIndex].pop_front();

    SampleIndexBuffer[bufferIndex].shrink_to_fit();
    LOPhaseBuffer[bufferIndex].shrink_to_fit();
    PWFreqBuffer[bufferIndex].shrink_to_fit();
    PWPhaseBuffer[bufferIndex].shrink_to_fit();
    PWMagBuffer[bufferIndex].shrink_to_fit();
    PatchVoltageBuffer[bufferIndex].shrink_to_fit();

    
  }
  
  void PlaneWaveSignalGenerator::InitializeBuffers(unsigned fieldbuffersize)
  {

    const unsigned nchannels = fNChannels;
    const int nReceivers = fNPatchesPerStrip;

    
    FieldBuffer aFieldBuffer;

    SampleIndexBuffer = aFieldBuffer.InitializeUnsignedBuffer(nchannels, nReceivers, fieldbuffersize);
    LOPhaseBuffer = aFieldBuffer.InitializeBuffer(nchannels, nReceivers, fieldbuffersize);
    PWFreqBuffer = aFieldBuffer.InitializeBuffer(nchannels, nReceivers, fieldbuffersize);
    PWMagBuffer = aFieldBuffer.InitializeBuffer(nchannels, nReceivers, fieldbuffersize);
    PWPhaseBuffer = aFieldBuffer.InitializeBuffer(nchannels, nReceivers, fieldbuffersize);   
    PatchVoltageBuffer = aFieldBuffer.InitializeBuffer(nchannels, nReceivers, fieldbuffersize);
    
    
  }

  void PlaneWaveSignalGenerator::InitializePatchArray()
  {

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

	    modelPatch.SetCenterPosition({patchRadius * cos(theta) , patchRadius * sin(theta) , zPosition }); 
	    modelPatch.SetPolarizationDirection({sin(theta), -cos(theta), 0.}); 
	    modelPatch.SetNormalDirection({-cos(theta), -sin(theta), 0.}); //Say normals point inwards
	    allChannels[channelIndex].AddReceiver(modelPatch);
	  }
      }
  }



  bool PlaneWaveSignalGenerator::DoGenerate( Signal* aSignal )
  {

    InitializePatchArray();

    // initialize FIR filter array
    for (unsigned i=0; i < sizeof(FIR_array)/sizeof(FIR_array[0]); i++)
      {  
	FIR_array[i] = 0.;
      }
    
    ProcessFIRFilter(1);
    nfilterbins = GetNFilterBins();
    int tempfieldbuffersize = 100;
 
    InitializeBuffers(tempfieldbuffersize);

    //text file for VI for testing.
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
    //    voltagefile.close();


    return true;
  }

} /* namespace locust */

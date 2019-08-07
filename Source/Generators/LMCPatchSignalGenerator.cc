/*
 * LMCPatchSignalGenerator.cc
 *
 *  Created on: Feb 19, 2018
 *      Author: pslocum
 */

#include "LMCPatchSignalGenerator.hh"
#include "LMCEventHold.hh"
#include "LMCRunKassiopeia.hh"

#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>

#include "LMCGlobalsDeclaration.hh"
#include "LMCDigitizer.hh"
#include <chrono>


namespace locust
{
    LOGGER( lmclog, "PatchSignalGenerator" );

    MT_REGISTER_GENERATOR(PatchSignalGenerator, "patch-signal");

    PatchSignalGenerator::PatchSignalGenerator( const std::string& aName ) :
        Generator( aName ),
        fLO_Frequency( 0.),
        fArrayRadius( 0. ),
        fNPatchesPerStrip( 0. ),
        fPatchSpacing( 0. ),
        fPowerCombiner( 0 ),
		fRJunction( 0.3 ),
        gxml_filename("blank.xml"),
		fTextFileWriting( 0 ),
        phiLO_t(0.),
        VoltagePhase_t {0.},
		gfilter_filename("blank.txt"),
        fFilter_resolution( 0. ),
        EFieldBuffer( 1 ),
        EPhaseBuffer( 1 ),
        EAmplitudeBuffer( 1 ),
        EFrequencyBuffer( 1 ),
        LOPhaseBuffer( 1 ),
        IndexBuffer( 1 ),
        PatchFIRBuffer( 1 ),
        fFieldBufferSize( 50 ),
        fFieldBufferMargin( 25 )

    {
        fRequiredSignalState = Signal::kTime;
    }

    PatchSignalGenerator::~PatchSignalGenerator()
    {
    }

    bool PatchSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam->has( "filter-filename" ) )
        {
            gfilter_filename = aParam->get_value< std::string >( "filter-filename" );
        }

        if( aParam->has( "filter-resolution" ) )
        {
            fFilter_resolution = aParam->get_value< double >( "filter-resolution" );
        }

        if( aParam->has( "buffer-size" ) )
        {
        	fFieldBufferSize = aParam->get_value< double >( "buffer-size" );
        }

        if( aParam->has( "buffer-margin" ) )
        {
        	fFieldBufferMargin = aParam->get_value< double >( "buffer-margin" );
        }



        if( aParam == NULL) return true;

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
        if( aParam->has( "xml-filename" ) )
        {
            gxml_filename = aParam->get_value< std::string >( "xml-filename" );
        }
        if( aParam->has( "text-filewriting" ) )
        {
            fTextFileWriting = aParam->get_value< bool >( "text-filewriting" );
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
    	else if (aParam->get_value< std::string >( "feed") == "voltage-divider")
    	  fPowerCombiner = 5;
    	else
    	  fPowerCombiner = 0;  // default
          }
        if( aParam->has( "junction-resistance" ) )
        {
            fRJunction = aParam->get_value< double >( "junction-resistance" );
        }


        return true;
    }

    void PatchSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    static void* KassiopeiaInit(const std::string &aFile)
    {
        //RunKassiopeia *RunKassiopeia1 = new RunKassiopeia;
        RunKassiopeia RunKassiopeia1;
        RunKassiopeia1.Run(aFile);
        RunKassiopeia1.~RunKassiopeia();
        //delete RunKassiopeia1;

        return 0;
    }



    static void WakeBeforeEvent()
    {
        fPreEventCondition.notify_one();
        return;
    }

    static bool ReceivedKassReady()
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        printf("LMC about to wait ..\n");

        if( !fKassEventReady)
        {
            std::unique_lock< std::mutex >tLock( fKassReadyMutex );
            fKassReadyCondition.wait( tLock );
        }

        if (fFalseStartKassiopeia)  // workaround for some Macs
        {
            std::unique_lock< std::mutex >tLock( fKassReadyMutex );
            fKassReadyCondition.wait( tLock );
        }

        return true;
    }



    double PatchSignalGenerator::GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition )
    {
        return aReceiverTime - aParticleTime - (aReceiverPosition - aParticlePosition).Magnitude() / LMCConst::C();
    }


    double GetPatchStepRoot(const locust::Particle aParticle, double aReceiverTime, LMCThreeVector aReceiverPosition, double aSpaceTimeInterval)
    {
        double tRetardedTime = aParticle.GetTime(true);
        return tRetardedTime + aSpaceTimeInterval;
    }


    double GetMismatchFactor(double f)  
    {
        //f /= 2.*LMCConst::Pi();
        // placeholder = 1 - mag(S11)
        // fit to HFSS output
        //double MismatchFactor = 1. - (-5.39e16 / ((f-25.9141e9)*(f-25.9141e9) + 7.23e16) + 0.88);
        //    printf("dopplerfrequency is %f and mismatchfactor is %g\n", f, MismatchFactor);  getchar();
        double MismatchFactor = 0.85;  // punt.
        return MismatchFactor;
    }

    double GetAOIFactor(LMCThreeVector IncidentKVector, double PatchPhi)
    {
        LMCThreeVector PatchNormalVector;
        PatchNormalVector.SetComponents(cos(PatchPhi), sin(PatchPhi), 0.0);
        double AOIFactor = fabs(IncidentKVector.Unit().Dot(PatchNormalVector));
        //printf("cos aoi is %f\n", AOIFactor);
        return AOIFactor;
    }


    // fields incident on patch.
    void RecordIncidentFields(FILE *fp, LMCThreeVector IncidentMagneticField, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double PatchPhi, double DopplerFrequency)
    {
        double AOIFactor = GetAOIFactor(IncidentKVector, PatchPhi);  // k dot patchnormal
        LMCThreeVector PatchPolarizationVector;
        PatchPolarizationVector.SetComponents(-sin(PatchPhi), cos(PatchPhi), 0.0);
        LMCThreeVector PatchCrossPolarizationVector;
        PatchCrossPolarizationVector.SetComponents(0.0, 0.0, 1.0);  // axial direction.

        double EFieldPatch = IncidentElectricField.Dot(PatchPolarizationVector);
        double BFieldPatch = IncidentMagneticField.Dot(PatchCrossPolarizationVector);
        fprintf(fp, "%10.4g %10.4g %10.4g %10.4g\n", EFieldPatch, BFieldPatch, DopplerFrequency/2./LMCConst::Pi(), t_old);
    }


    double* PatchSignalGenerator::GetFIRFilter(int nskips)
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

    int PatchSignalGenerator::GetNFilterBins(double* filterarray)
    {
    int nbins = 0;
    for (int i=0; i<1000; i++)
      {
      if (filterarray[i]>0.) nbins += 1;
      }
    return nbins;
    }


    double PatchSignalGenerator::GetFIRSample(double* filterarray, int nfilterbins, double dtfilter, unsigned channel, unsigned patch, double AcquisitionRate)
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

   	for (int i=0; i < nfilterbins; i++)  // populate filter with field.
      {
    	  HilbertPhase += 2.*3.1415926*fieldfrequency*dtfilter;
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

//    return EFieldBuffer[channel*fNPatchesPerStrip+patch].front(); // debug
//    return HilbertMag*cos(HilbertPhase);  // debug

    }



    // EField cross pol with aoi dot product, at patch.
    double GetEFieldCoPol(PatchAntenna* currentPatch, LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double PatchPhi, double DopplerFrequency)
    {
        double AOIFactor = GetAOIFactor(IncidentKVector, PatchPhi);  // k dot patchnormal
        LMCThreeVector PatchPolarizationVector = currentPatch->GetPolarizationDirection();
        double EFieldCoPol = IncidentElectricField.Dot(PatchPolarizationVector) * AOIFactor;

        return EFieldCoPol;
    }


    // voltage amplitude induced at patch.
    double GetVoltageAmplitude(LMCThreeVector IncidentElectricField, LMCThreeVector IncidentKVector, double PatchPhi, double DopplerFrequency)
    {
        double AntennaFactor = 1./400.;
        double MismatchFactor = GetMismatchFactor(DopplerFrequency);
        double AOIFactor = GetAOIFactor(IncidentKVector, PatchPhi);  // k dot patchnormal
        LMCThreeVector PatchPolarizationVector;
        PatchPolarizationVector.SetComponents(-sin(PatchPhi), cos(PatchPhi), 0.0);
        double VoltageAmplitude = fabs( AntennaFactor * IncidentElectricField.Dot(PatchPolarizationVector) * MismatchFactor * AOIFactor);
        //double VoltageAmplitude = fabs( AntennaFactor * IncidentElectricField.Magnitude()); // test case.  

        //    if (VoltageAmplitude>0.) {printf("IncidentElectricField.Dot(PatchPolarizationVector) is %g and VoltageAmplitude is %g\n", IncidentElectricField.Dot(PatchPolarizationVector), VoltageAmplitude); getchar();}
        return VoltageAmplitude;
    }


    void PatchSignalGenerator::AddOneFIRVoltageToStripSum(Signal* aSignal, double VoltageFIRSample, double phi_LO, unsigned channelIndex, unsigned patchIndex)
    {

    	PowerCombiner aPowerCombiner;

        if (fPowerCombiner == 0 ) //corporate feed, for testing
          {
    	    VoltageFIRSample *= aPowerCombiner.GetCorporateVoltageDamping();
          }

        if (fPowerCombiner == 3) // seven-eighths power combining, center fed strip
          {
    	    VoltageFIRSample *= aPowerCombiner.GetSevenEighthsVoltageDamping(fNPatchesPerStrip, patchIndex);
          }
        if (fPowerCombiner == 5)
          {
            VoltageFIRSample *= aPowerCombiner.GetVoltageDividerWeight(fRJunction, 1.0, 10.e6, fNPatchesPerStrip, patchIndex);
          }

    	aSignal->LongSignalTimeComplex()[IndexBuffer[channelIndex*fNPatchesPerStrip+patchIndex].front()][0] += 2.*VoltageFIRSample * sin(phi_LO);
        aSignal->LongSignalTimeComplex()[IndexBuffer[channelIndex*fNPatchesPerStrip+patchIndex].front()][1] += 2.*VoltageFIRSample * cos(phi_LO);


    }


    // z-index ranges from 0 to npatches-per-strip-1.
    void PatchSignalGenerator::AddOnePatchVoltageToStripSum(Signal* aSignal, double VoltageAmplitude, double VoltagePhase, double phi_LO, unsigned channelindex, unsigned z_index, double DopplerFrequency)
    {
    	PowerCombiner aPowerCombiner;
        if (fPowerCombiner == 1)  // series feed
        {
        	//lossless series feed with amp at one end:
        	//VoltagePhase += aPowerCombiner.GetLinePhaseCorr(z_index, DopplerFrequency);
        }
        if (fPowerCombiner == 2) // quadrature feed
        {
        	// assume 2PI delay between junctions, so we don't calculated phase mismatches.
        	// instead calculate damping on voltage amplitude:
            int njunctions = (int)fabs(z_index - fNPatchesPerStrip/2);
           // VoltageAmplitude *= aPowerCombiner.GetVoltageDamping(fNPatchesPerStrip, z_index);
        }

//	        if (VoltageAmplitude>0.) {printf("voltageamplitude is %g\n", VoltageAmplitude); getchar();}
        aSignal->LongSignalTimeComplex()[channelindex][0] += VoltageAmplitude * cos(VoltagePhase - phi_LO);
        aSignal->LongSignalTimeComplex()[channelindex][1] += VoltageAmplitude * sin(VoltagePhase - phi_LO);
	//        if (VoltageAmplitude>0.) {printf("summedvoltageamplitude is %g\n", aSignal->LongSignalTimeComplex()[channelindex][0]); getchar();}                           


    }


    void* PatchSignalGenerator::DriveAntenna(FILE *fp, int PreEventCounter, unsigned index, Signal* aSignal, double* filterarray, unsigned nfilterbins, double dtfilter)
    {
        if (PreEventCounter > 0)  // new event starting.                                                    
        {
            // initialize patch voltage phases.                                                             
            for (unsigned i=0; i < sizeof(VoltagePhase_t)/sizeof(VoltagePhase_t[0]); i++)
            {
                VoltagePhase_t[i] = {0.};
            }
        }

        locust::Particle tCurrentParticle = fParticleHistory.back();
        int CurrentIndex;
        const int signalSize = aSignal->TimeSize();

        const double kassiopeiaTimeStep = fabs(fParticleHistory[0].GetTime() - fParticleHistory[1].GetTime());
        const int historySize = fParticleHistory.size();
        unsigned sampleIndex = 0;

        //Receiver Properties
        phiLO_t += 2. * LMCConst::Pi() * fLO_Frequency * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
        double tReceiverTime = t_old;
        double tRetardedTime = 0.; //Retarded time of particle corresponding to when emission occurs, reaching receiver at tReceiverTime

        double tSpaceTimeInterval=99.;
        double dtRetarded=0;
        double tTolerance=1e-23;

        PatchAntenna *currentPatch;

        for(int channelIndex = 0; channelIndex < allChannels.size(); ++channelIndex)
        {
            double PatchPhi = (double)channelIndex*360./allChannels.size()*LMCConst::Pi()/180.; // radians.    
            for(int patchIndex = 0; patchIndex < allChannels[channelIndex].size(); ++patchIndex)
            {
                currentPatch = &allChannels[channelIndex][patchIndex]; 
                sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample

                if(fParticleHistory.front().GetTime()<=3.*kassiopeiaTimeStep)
                {
                    fParticleHistory.front().Interpolate(0);
                    if(GetSpaceTimeInterval(fParticleHistory.front().GetTime(true), tReceiverTime , fParticleHistory.front().GetPosition(true), currentPatch->GetPosition() ) < 0 )
                    {
                        //printf("Skipping! out of Bounds!: tReceiverTime=%e\n",tReceiverTime);
                        continue;
                    }
                }

                if(currentPatch->GetPreviousRetardedIndex() == -99.)
                {
                    CurrentIndex=FindNode(tReceiverTime);
                    tCurrentParticle = fParticleHistory[CurrentIndex];
                    tRetardedTime = tReceiverTime - (tCurrentParticle.GetPosition() - currentPatch->GetPosition() ).Magnitude() /  LMCConst::C();
                }
                else
                {
                    CurrentIndex = currentPatch->GetPreviousRetardedIndex();
                    tRetardedTime = currentPatch->GetPreviousRetardedTime() + 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
                }

                CurrentIndex = FindNode(tRetardedTime);
                CurrentIndex = std::min(std::max(CurrentIndex,0) , historySize - 1);

                tCurrentParticle = fParticleHistory[CurrentIndex];
                tCurrentParticle.Interpolate(tRetardedTime);
                tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tReceiverTime, tCurrentParticle.GetPosition(true), currentPatch->GetPosition());

                double tOldSpaceTimeInterval=99.;

                //Converge to root
                for(int j=0;j<25;++j)
                {
                    tRetardedTime = GetPatchStepRoot(tCurrentParticle, tReceiverTime, currentPatch->GetPosition(), tSpaceTimeInterval);
                    tCurrentParticle.Interpolate(tRetardedTime);

                    //Change the kassiopeia step we expand around if the interpolation time displacement is too large
                    if(fabs(tCurrentParticle.GetTime(true) - tCurrentParticle.GetTime(false)) > kassiopeiaTimeStep)
                    {
                        CurrentIndex=FindNode(tRetardedTime);
                        tCurrentParticle=fParticleHistory[CurrentIndex];
                        tCurrentParticle.Interpolate(tRetardedTime);
                    }

                    tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tReceiverTime, tCurrentParticle.GetPosition(true), currentPatch->GetPosition());
                    tOldSpaceTimeInterval = tSpaceTimeInterval;
                }


                currentPatch->SetPreviousRetardedIndex(CurrentIndex);
                currentPatch->SetPreviousRetardedTime(tRetardedTime);

                LMCThreeVector tDirection = currentPatch->GetPosition() - tCurrentParticle.GetPosition(true);
                double tVelZ = tCurrentParticle.GetVelocity(true).Z();
                double tCosTheta =  tVelZ * tDirection.Z() /  tDirection.Magnitude() / fabs(tVelZ);
                double tDopplerFrequency  = tCurrentParticle.GetCyclotronFrequency() / ( 1. - fabs(tVelZ) / LMCConst::C() * tCosTheta);


                if (VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex]>0.)  // not first sample                                                        
                {
                    VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex] += tDopplerFrequency * 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());
                }
                else  // if this is the first light at this patch, the voltage phase doesn't advance for the full dt.
                {              
                    VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex] += tDopplerFrequency * tRetardedTime;
                    //printf("tDopplerFrequency is %g\n", tDopplerFrequency); getchar();
                }      

 		        double tEFieldCoPol = GetEFieldCoPol(currentPatch, tCurrentParticle.CalculateElectricField(currentPatch->GetPosition()), tCurrentParticle.CalculateElectricField(currentPatch->GetPosition()).Cross(tCurrentParticle.CalculateMagneticField(currentPatch->GetPosition())), PatchPhi, tDopplerFrequency);
                if (fTextFileWriting==1) RecordIncidentFields(fp, tCurrentParticle.CalculateMagneticField(currentPatch->GetPosition()), tCurrentParticle.CalculateElectricField(currentPatch->GetPosition()), tCurrentParticle.CalculateElectricField(currentPatch->GetPosition()).Cross(tCurrentParticle.CalculateMagneticField(currentPatch->GetPosition())) , PatchPhi, tDopplerFrequency);

 	            FillBuffers(aSignal, tDopplerFrequency, tEFieldCoPol, phiLO_t, index, channelIndex, patchIndex, 0);
 	            double VoltageFIRSample = GetFIRSample(filterarray, nfilterbins, dtfilter, channelIndex, patchIndex, fAcquisitionRate*aSignal->DecimationFactor());
 	            AddOneFIRVoltageToStripSum(aSignal, VoltageFIRSample, phiLO_t, channelIndex, patchIndex);
                PopBuffers(channelIndex, patchIndex);

// 	            AddOnePatchVoltageToStripSum(aSignal, tVoltageAmplitude, VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex], phiLO_t, sampleIndex, patchIndex, tDopplerFrequency);

            } // patch loop

        } // channels loop

        t_old += 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

        return 0;
    }

    //Return index of fParticleHistory particle closest to the time we are evaluating
    int PatchSignalGenerator::FindNode(double tNew) const
    {
        std::deque<locust::Particle>::iterator it;
        it = std::upper_bound( fParticleHistory.begin() , fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );

        int tNodeIndex = it - fParticleHistory.begin();

        return tNodeIndex;
    }



    void PatchSignalGenerator::FillBuffers(Signal* aSignal, double DopplerFrequency, double EFieldValue, double LOPhase, unsigned index, unsigned channel, unsigned patch, unsigned dtauConvolutionTime)
    {
    EFieldBuffer[channel*fNPatchesPerStrip+patch].push_back(EFieldValue);
    EFrequencyBuffer[channel*fNPatchesPerStrip+patch].push_back(DopplerFrequency/2./LMCConst::Pi());
    LOPhaseBuffer[channel*fNPatchesPerStrip+patch].push_back(LOPhase);
    IndexBuffer[channel*fNPatchesPerStrip+patch].push_back(channel*aSignal->TimeSize()*aSignal->DecimationFactor() + index);

    }





    void PatchSignalGenerator::PopBuffers(unsigned channel, unsigned patch)
    {

    	EFieldBuffer[channel*fNPatchesPerStrip+patch].pop_front();
    	EFrequencyBuffer[channel*fNPatchesPerStrip+patch].pop_front();
    	LOPhaseBuffer[channel*fNPatchesPerStrip+patch].pop_front();
    	IndexBuffer[channel*fNPatchesPerStrip+patch].pop_front();

    	EFieldBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        EFrequencyBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        LOPhaseBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();
        IndexBuffer[channel*fNPatchesPerStrip+patch].shrink_to_fit();

    }




    void PatchSignalGenerator::InitializeBuffers(unsigned filterbuffersize, unsigned fieldbuffersize)
    {

    FieldBuffer aFieldBuffer;

    EFieldBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
    EFrequencyBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
    LOPhaseBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);
    IndexBuffer = aFieldBuffer.InitializeUnsignedBuffer(fNChannels, fNPatchesPerStrip, fieldbuffersize);

    PatchFIRBuffer = aFieldBuffer.InitializeBuffer(fNChannels, fNPatchesPerStrip, filterbuffersize);


    }


    void PatchSignalGenerator::CleanupBuffers()
    {
    FieldBuffer aFieldBuffer;
    EFieldBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    EFrequencyBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    LOPhaseBuffer = aFieldBuffer.CleanupBuffer(EFieldBuffer);
    IndexBuffer = aFieldBuffer.CleanupBuffer(IndexBuffer);

    }


  double RotateZ(int component, double angle, double x, double y)
    {
      double newcomponent = 0.;
      if (component==0)
        {
        newcomponent = x*cos(angle) - y*sin(angle);
        }
      if (component==1)
        {
        newcomponent = x*sin(angle) + y*cos(angle);
        }

      return newcomponent;
    }


    void PatchSignalGenerator::InitializePatchArray()
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



    bool PatchSignalGenerator::DoGenerate( Signal* aSignal )
    {

        FILE *fp = fopen("incidentfields.txt", "w");

        InitializePatchArray();


        //n samples for event spacing.
        int PreEventCounter = 0;
        const int NPreEventSamples = 150000;
        fKassTimeStep = 1./(fAcquisitionRate*1.e6*aSignal->DecimationFactor());

        std::thread Kassiopeia(KassiopeiaInit, gxml_filename);     // spawn new thread
        fRunInProgress = true;

        double* filterarray = GetFIRFilter(1);
        unsigned nfilterbins = GetNFilterBins(filterarray);
        unsigned nfieldbufferbins = fFieldBufferSize;
        double dtfilter = fFilter_resolution;
        unsigned dtauConvolutionTime = 0;
        InitializeBuffers(nfilterbins, nfieldbufferbins);

        for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        {
            if ((!fEventInProgress) && (fRunInProgress) && (!fPreEventInProgress))
            {
                if (ReceivedKassReady()) fPreEventInProgress = true;
            }

            if (fPreEventInProgress)
            {
                PreEventCounter += 1;
                //printf("preeventcounter is %d\n", PreEventCounter);
                if (PreEventCounter > NPreEventSamples)  // finished noise samples.  Start event.
                {
                    fPreEventInProgress = false;  // reset.
                    fEventInProgress = true;
                    //printf("LMC about to wakebeforeevent\n");
                    WakeBeforeEvent();  // trigger Kass event.
                }
            }

            if (fEventInProgress)  // fEventInProgress
                if (fEventInProgress)  // check again.
                {
                    //printf("waiting for digitizer trigger ... index is %d\n", index);
                    std::unique_lock< std::mutex >tLock( fMutexDigitizer, std::defer_lock );
                    tLock.lock();
                    fDigitizerCondition.wait( tLock );
                    if (fEventInProgress)
                    {
                        //printf("about to drive antenna, PEV is %d\n", PreEventCounter);
                        DriveAntenna(fp, PreEventCounter, index, aSignal, filterarray, nfilterbins, dtfilter);
                        PreEventCounter = 0; // reset
                    }
                    tLock.unlock();
                }

        }  // for loop

        printf("finished signal loop\n");

        fclose(fp);
        CleanupBuffers();
        fRunInProgress = false;  // tell Kassiopeia to finish.
        fDoneWithSignalGeneration = true;  // tell LMCCyclotronRadExtractor
        //if (fEventInProgress)
        //  if (ReceivedKassReady())
        WakeBeforeEvent();
        Kassiopeia.join();

        return true;
    }

} /* namespace locust */


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
        fNPatchesPerStrip( 0. ),
        fPatchSpacing( 0. ),
        fPowerCombiner( 0 ),
        phiLO_t(0.),
        VoltagePhase_t {0.}
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
            else if (aParam->get_value< std::string >( "feed" ) == "quadrature")
                fPowerCombiner = 2;
            else
            	fPowerCombiner = 0;  // default
        }

        return true;
    }

    void PlaneWaveSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }



    double PlaneWaveSignalGenerator::GetMismatchFactor(double f)
    {
        //f /= 2.*LMCConst::Pi();
        // placeholder = 1 - mag(S11)
        // fit to HFSS output
        //double MismatchFactor = 1. - (-5.39e16 / ((f-25.9141e9)*(f-25.9141e9) + 7.23e16) + 0.88);
        //    printf("dopplerfrequency is %f and mismatchfactor is %g\n", f, MismatchFactor);  getchar();
        double MismatchFactor = 0.85;  // punt.
        return MismatchFactor;
    }

    double PlaneWaveSignalGenerator::GetAOIFactor(LMCThreeVector IncidentKVector, double PatchPhi) const
    {
        LMCThreeVector PatchNormalVector;
        PatchNormalVector.SetComponents(cos(PatchPhi), sin(PatchPhi), 0.0);
        double AOIFactor = fabs(IncidentKVector.Unit().Dot(PatchNormalVector));
        //printf("cos aoi is %f\n", AOIFactor);
        return AOIFactor;
    }

    double PlaneWaveSignalGenerator::GetVoltageAmpFromPlaneWave()
    {
        double AntennaFactor = 1./400.;

        // S = epsilon0 c E0^2 / 2.  // power/area
        //  0.6e-21 W/Hz * 24.e3 Hz / (0.00375*0.002916) = S = 1.3e-12 W/m^2
        // We should detect 0.6e-21 W/Hz * 24.e3 Hz in Katydid.
        // E0 = sqrt(2.*S/epsilon0/c)
        // effective patch area 0.00004583662 m^2 

        double S = 0.6e-21*24.e3/(0.00004271);  // W/m^2, effective aperture.
        double E0 = sqrt(2.*S/(LMCConst::EpsNull() * LMCConst::C()));
        //    double E0 = 1.0; // V/m, test case
        double amplitude = E0*AntennaFactor;  // volts
        return amplitude;


    }



    // z-index ranges from 0 to npatches-per-strip-1.
    void PlaneWaveSignalGenerator::AddOnePatchVoltageToStripSum(Signal* aSignal, double VoltageAmplitude, double VoltagePhase, double phi_LO, unsigned channelindex, unsigned z_index, double DopplerFrequency)
    {
    	PowerCombiner aPowerCombiner;
        if (fPowerCombiner == 1)  // series feed
        {
        	//lossless series feed with amp at one end:
        	VoltagePhase += aPowerCombiner.GetLinePhaseCorr(z_index, DopplerFrequency);
        }
        if (fPowerCombiner == 2) // quadrature feed
        {
            if (fPhaseDelay) // parameter from json
            {
	      VoltagePhase += aPowerCombiner.GetCenterFedLinePhaseCorr(fNPatchesPerStrip, z_index, DopplerFrequency, fPatchSpacing);
            }
        	// instead calculate damping on voltage amplitude:
            int njunctions = (int)fabs(z_index - fNPatchesPerStrip/2);
            VoltageAmplitude *= aPowerCombiner.GetVoltageDamping(njunctions);
        }
        aSignal->LongSignalTimeComplex()[channelindex][0] += VoltageAmplitude * cos(VoltagePhase - phi_LO);
        aSignal->LongSignalTimeComplex()[channelindex][1] += VoltageAmplitude * sin(VoltagePhase - phi_LO);
//	        if (VoltageAmplitude>0.) {printf("summedvoltageamplitude is %g\n", aSignal->LongSignalTimeComplex()[channelindex][0]); getchar();}                           


    }


    void* PlaneWaveSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal)
    {

        const int signalSize = aSignal->TimeSize();
        unsigned sampleIndex = 0;
        double testphase = 0.;

        //Receiver Properties
        phiLO_t += 2. * LMCConst::Pi() * fLO_Frequency / (1.e6*fAcquisitionRate*aSignal->DecimationFactor());

        PatchAntenna *currentPatch;

        for(int channelIndex = 0; channelIndex < allChannels.size(); ++channelIndex)
        {
            for(int patchIndex = 0; patchIndex < allChannels[channelIndex].size(); ++patchIndex)
            {
                currentPatch = &allChannels[channelIndex][patchIndex]; 
                sampleIndex = channelIndex*signalSize*aSignal->DecimationFactor() + index;  // which channel and which sample

                VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex] +=
 2. * LMCConst::Pi() * fRF_Frequency / (1.e6 * fAcquisitionRate * aSignal->DecimationFactor()); // phi =+ f*dt
//            printf("voltagephase_t is %g\n", VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex]); getchar();

                double tVoltageAmplitude = GetVoltageAmpFromPlaneWave();
                AddOnePatchVoltageToStripSum(aSignal, tVoltageAmplitude, VoltagePhase_t[channelIndex*fNPatchesPerStrip+patchIndex], phiLO_t, sampleIndex, patchIndex, fLO_Frequency);

            } // patch loop

        } // channels loop

        return 0;
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


        //n samples for event spacing.
        int PreEventCounter = 0;
        const int NPreEventSamples = 150000;
        PreEventCounter = NPreEventSamples; // jump past wait time.
        for (unsigned i=0; i < sizeof(VoltagePhase_t)/sizeof(VoltagePhase_t[0]); i++)
            {  // initialize voltage phases.
                VoltagePhase_t[i] = {0.0};
            }



        for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
          {
          DriveAntenna(PreEventCounter, index, aSignal);
          }

        return true;
    }

} /* namespace locust */

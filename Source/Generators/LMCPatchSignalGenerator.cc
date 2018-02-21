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
#include "LMCHFSSReader.hh"
#include "LMCSimulationController.hh"
#include <chrono>


namespace locust
{
    LOGGER( lmclog, "PatchSignalGenerator" );

    MT_REGISTER_GENERATOR(PatchSignalGenerator, "patch-signal");

    PatchSignalGenerator::PatchSignalGenerator( const std::string& aName ) :
            Generator( aName ),
            fWriteNFD(0.),
            fLO_Frequency( 0.),
            gxml_filename("blank.xml")
    {
        fRequiredSignalState = Signal::kTime;
    }

    PatchSignalGenerator::~PatchSignalGenerator()
    {
    }

    bool PatchSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;

        if( aParam->has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam->get_value< double >( "lo-frequency" );
        }
        if( aParam->has( "xml-filename" ) )
        {
            gxml_filename = aParam->get_value< std::string >( "xml-filename" );
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
    	RunKassiopeia *RunKassiopeia1 = new RunKassiopeia;
    	RunKassiopeia1->Run(aFile);
    	delete RunKassiopeia1;

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
            printf("LMC Got the fKassReadyCondition signal\n");
        }

        return true;
    }



    double PatchSignalGenerator::GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const LMCThreeVector &aParticlePosition, const LMCThreeVector &aReceiverPosition )
    {
        //return pow(aReceiverTime - aParticleTime,2.) - (aReceiverPosition - aParticlePosition).MagnitudeSquared() / pow(LMCConst::C() , 2.);
        return aReceiverTime - aParticleTime - (aReceiverPosition - aParticlePosition).Magnitude() / LMCConst::C();
    }

    double PatchGetStepRoot(const locust::Particle aParticle, double aReceiverTime, LMCThreeVector aReceiverPosition, double aSpaceTimeInterval, const int aStepOrder = 0)
    {
        double tRetardedTime = aParticle.GetTime(true); //interpolate!!!

        double c=LMCConst::C();

        if(aStepOrder==0)
        {
            double tCorrection = sqrt(fabs(aSpaceTimeInterval));
            double tSign;
            (aSpaceTimeInterval > 0) ? tSign = 1. : tSign = -1.;
            
            return tRetardedTime + aSpaceTimeInterval;
            //return tRetardedTime + tSign * aSpaceTimeInterval;
        }

        LMCThreeVector tNewPosition = aParticle.GetPosition(true);
        LMCThreeVector tNewVelocity = aParticle.GetVelocity(true);

        LMCThreeVector tReceiverVector = aReceiverPosition - tNewPosition;
        double tReceiverDistance = tReceiverVector.Magnitude();

        //Newtons Method X_{n+1} = X_{n} - f(X_{n}) / f'(X_{n})
        double fZero=pow((aReceiverTime - tRetardedTime),2.)-pow(tReceiverDistance,2.)/(c*c);
        double fZeroPrime=2.*((tRetardedTime-aReceiverTime)-tNewVelocity.Dot(tNewPosition)/(c*c)+tNewVelocity.Dot(aReceiverPosition)/(c*c));
        double tNewtonRatio = fZero / fZeroPrime;

        if(aStepOrder==1)
        {
            return tRetardedTime-tNewtonRatio;
        }

        //Householders Method
        LMCThreeVector tNewAcceleration = aParticle.GetAcceleration(true);
        double fZeroDoublePrime = 2. * (1. - tNewVelocity.Dot(tNewVelocity)/(c*c)-tNewAcceleration.Dot(tNewPosition-aReceiverPosition)/(c*c));

        if(aStepOrder==2)
        {
            return tRetardedTime-tNewtonRatio* ( 1. + ( tNewtonRatio * fZeroDoublePrime) / ( 2. * fZeroPrime));
        }

        LERROR( lmclog, "Need to put root finding method with order 0-2!" );
        return 0;
    }


    double GetVoltageAmplitude()
    {
    	double VoltageAmplitude = 1.e-8; // placeholder
    	return VoltageAmplitude;
    }


    void AddOnePatchVoltageToStripSum(Signal* aSignal, double VoltageAmplitude, double VoltagePhase, double phi_LO, unsigned channelindex)
    {
        aSignal->LongSignalTimeComplex()[channelindex][0] += VoltageAmplitude * cos(VoltagePhase - phi_LO);
        aSignal->LongSignalTimeComplex()[channelindex][1] += VoltageAmplitude * sin(VoltagePhase - phi_LO);
    }



    void* PatchSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal)
    {

        locust::Particle tCurrentParticle = fParticleHistory.back();
        int CurrentIndex;
        HFSSReader HFRead;


    	SimulationController SimulationController1;
        const unsigned nchannels = SimulationController1.GetNChannels();
        const double dx=0.01;

        const int nelementsZ = 9;

        double theta = 0.;  // azimuthal angle of each antenna element.
        double radius = 0.05; // radius of each antenna element in xy plane.
        static double tVoltagePhase[10000] = {0.};  // this is not resetting at the beginning of each event.  big problem.
        static double phi_LO = 0.;

        const int signalSize = aSignal->TimeSize();
        unsigned patchindex = 0;
        unsigned channelindex = 0;

        //Receiver Properties
        double tReceiverTime = t_old;
        LMCThreeVector tReceiverPosition;

        double tRetardedTime = 0.; //Retarded time of particle corresponding to when emission occurs, reaching receiver at tReceiverTime
        double tTotalPower=0.;

        double tSpaceTimeInterval=99.;
        double dtRetarded=0;
        double tTolerance=1e-23;
        const double dtStepSize = fabs(fParticleHistory[0].GetTime() - fParticleHistory[1].GetTime());

        const int HistorySize = fParticleHistory.size();

        phi_LO+= 2. * LMCConst::Pi() * fLO_Frequency * fDigitizerTimeStep;  // this has to happen outside the signal generating loop.

        //int tAverageIterations=0; //Performance tracker. Count number of iterations to converge....

        for (int ch=0; ch<nchannels; ch++)  // number of patch strips and amplifier channels.
        {
        for (int z_index = 4; z_index<5; z_index++) // step through patch elements along z.  fix this.
        {
        // position patches in space:
        rReceiver = HFRead.GeneratePlane({dx,dx},7);//Argumemts: Size, resolution
        theta = (double)ch*360./nchannels*LMCConst::Pi()/180.;
        rReceiver = HFRead.RotateShift(rReceiver,{cos(theta),sin(theta),0.},{radius*cos(theta),radius*sin(theta),(double)(z_index-4)*0.01});//Arguments Normal vector, Position (m)
        PreviousTimes = std::vector<std::pair<int,double> >(rReceiver.size(),{-99.,-99.}); // initialize
        tTotalPower = 0.; // initialize
        patchindex = ch*nelementsZ + z_index;  // which patch element.
        channelindex = ch*signalSize*aSignal->DecimationFactor() + index;  // tells us which multichannel digitizer sample.


        for(unsigned i=0;i<rReceiver.size();++i)
        {
            tReceiverPosition = rReceiver[i];

            //Check if there is time for photon to reach receiver if particle is recently created
            if(fParticleHistory.front().GetTime()<=3.*dtStepSize)
            {
                fParticleHistory.front().Interpolate(0);
                if(GetSpaceTimeInterval(fParticleHistory.front().GetTime(true), tReceiverTime , fParticleHistory.front().GetPosition(true), tReceiverPosition) < 0 )
                {
                    //printf("Skipping! out of Bounds!: tReceiverTime=%e\n",tReceiverTime);
                    continue;
                }
            }

            if(PreviousTimes[i].first == -99.)
            {
                CurrentIndex=FindNode(tReceiverTime,dtStepSize,HistorySize-1);
                tCurrentParticle=fParticleHistory[CurrentIndex];

                tRetardedTime = tReceiverTime - (tCurrentParticle.GetPosition() - tReceiverPosition).Magnitude() /  LMCConst::C();
             }
            else
            {
                CurrentIndex = PreviousTimes[i].first;
                tRetardedTime = PreviousTimes[i].second + fDigitizerTimeStep;
            }

            CurrentIndex = std::min(std::max(CurrentIndex,0) , HistorySize - 1);

            CurrentIndex = FindNode(tRetardedTime,dtStepSize,CurrentIndex);

            tCurrentParticle = fParticleHistory[CurrentIndex];
            tCurrentParticle.Interpolate(tRetardedTime);
            tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tReceiverTime, tCurrentParticle.GetPosition(true), tReceiverPosition);

            double tOldSpaceTimeInterval=99.;


            //Converge to root
            for(int j=0;j<25;++j)
            {
                //++tAverageIterations;

                tRetardedTime = PatchGetStepRoot(tCurrentParticle, tReceiverTime, tReceiverPosition, tSpaceTimeInterval,0);

                tCurrentParticle.Interpolate(tRetardedTime);

                //Change the kassiopeia step we expand around if the interpolation time displacement is too large
                if(fabs(tCurrentParticle.GetTime(true) - tCurrentParticle.GetTime(false)) > dtStepSize)
                {
                    CurrentIndex=FindNode(tRetardedTime,dtStepSize,CurrentIndex);
                    tCurrentParticle=fParticleHistory[CurrentIndex];
                    tCurrentParticle.Interpolate(tRetardedTime);
                }
                //printf("%e %e\n",tCurrentParticle.GetTimeDisplacement(),dtStepSize*0.75);

                tSpaceTimeInterval = GetSpaceTimeInterval(tCurrentParticle.GetTime(true), tReceiverTime, tCurrentParticle.GetPosition(true), tReceiverPosition);

                tOldSpaceTimeInterval = tSpaceTimeInterval;
            }


            PreviousTimes[i].first = CurrentIndex;
            PreviousTimes[i].second = tRetardedTime;


            LMCThreeVector tECrossH = tCurrentParticle.CalculateElectricField(rReceiver[i]).Cross(tCurrentParticle.CalculateMagneticField(rReceiver[i]));
            LMCThreeVector tDirection = tReceiverPosition - tCurrentParticle.GetPosition(true);

            tTotalPower += dx * dx * tECrossH.Dot(tDirection.Unit()) / rReceiver.size() ;// * (fabs(tCurrentParticle.GetPosition(true).Z())<0.01);

            double tVelZ = tCurrentParticle.GetVelocity(true).Z();
            double tCosTheta =  tVelZ * tDirection.Z() /  tDirection.Magnitude() / fabs(tVelZ);
            double tDopplerFrequency  = tCurrentParticle.GetCyclotronFrequency() / ( 1. - fabs(tVelZ) / LMCConst::C() * tCosTheta);
            double tTransitTime = (tReceiverPosition - tCurrentParticle.GetPosition(true)).Magnitude() / LMCConst::C();

            if (tRetardedTime+tTransitTime > fDigitizerTimeStep)  // if the signal has been present for longer than fDigitizerTimeStep
              {
              tVoltagePhase[patchindex]+= tDopplerFrequency * fDigitizerTimeStep / rReceiver.size();
              }
            else  // if this is the first digitizer sample, the voltage phase doesn't advance for the full dt.
              {
                tVoltagePhase[patchindex]+=
                		tDopplerFrequency * (tReceiverTime - (tRetardedTime+tTransitTime)) / rReceiver.size();
//                printf("Retarding voltage phase:  fDigitizerTimeStep is %g and tReceiverTime is %g and tRetardedTime is %g\n and t_old is %g\n",
//                		fDigitizerTimeStep, tReceiverTime, tRetardedTime, t_old); getchar();
              }

            /*
            if (i==0)  // check receiver point 0 for each channel.  It should be the same each time.
            {
              printf("rx point 0:  zposition is %d and channel is %d, fcyc is %g and tVelZ is %g, dopplerfreq is %g, costheta is %f\n",
            		  z_index, ch, tDopplerFrequency, tCurrentParticle.GetCyclotronFrequency(), tVelZ, tCosTheta);
              printf("current index is %d\n", CurrentIndex);
            }
            */


        }  // i, rReceiver.size() loop.

        // By this point we will need to have the plane wave vector at the patch center.
        // We also now have the Doppler shifted frequency.
        // We need to calculate the voltage amplitude.
        double tVoltageAmplitude = GetVoltageAmplitude(/* parameters go here */);
        AddOnePatchVoltageToStripSum(aSignal, tVoltageAmplitude, tVoltagePhase[patchindex], phi_LO, channelindex);


        } // z_index waveguide element stepping loop.
        } // nchannels loop.

        t_old += fDigitizerTimeStep;

        return 0;
    }


    //Return iterator of fParticleHistory particle closest to the time we are evaluating
    int PatchSignalGenerator::FindNode(double tNew, double dtStepSize, int kIndexOld) const
    {
        int tHistorySize = fParticleHistory.size();

        //Make sure we are not out of bounds of array!!!
        kIndexOld = std::min( std::max(kIndexOld,0) , tHistorySize - 1 );

        double tOld = fParticleHistory[ kIndexOld ].GetTime();

        int kIndexMid=round((tNew-tOld)/dtStepSize) + kIndexOld;
        kIndexMid = std::min( std::max(kIndexMid,0) , tHistorySize - 1 );

        int kIndexSearchWidth;
        int kIndexRange[2];
        std::deque<locust::Particle>::iterator it;

        for(int i = 0 ; i < 15 ; ++i){

            kIndexSearchWidth = pow( 2 , i );
            kIndexRange[0] = kIndexMid - kIndexSearchWidth;
            kIndexRange[1] = kIndexMid + kIndexSearchWidth;

            kIndexRange[0] = std::max(kIndexRange[0], 0 );
            kIndexRange[1] = std::min(kIndexRange[1], tHistorySize - 1);

            if( tNew >= fParticleHistory[ kIndexRange[0] ].GetTime() && tNew <= fParticleHistory[ kIndexRange[1] ].GetTime())
            {
                //Get iterator pointing to particle step closest to tNew
                it = std::upper_bound( fParticleHistory.begin() , fParticleHistory.end() , tNew, [] (const double &a , const locust::Particle &b) { return a < b.GetTime();} );
                break;
            }
        }
        
        int tNodeIndex = it - fParticleHistory.begin();
        //if(tNodeIndex < 0 || tNodeIndex > (tHistorySize - 1))
        //{
        //    LERROR(lmclog, "OUT OF INDEX SEARCH");
        //}

        return tNodeIndex;
    }

    bool PatchSignalGenerator::DoGenerate( Signal* aSignal )
    {

    	//n samples for event spacing.
        int PreEventCounter = 0;
        const int NPreEventSamples = 150000;

        std::thread Kassiopeia (KassiopeiaInit, gxml_filename);     // spawn new thread
        fRunInProgress = true;

        for( unsigned index = 0; index < aSignal->DecimationFactor()*aSignal->TimeSize(); ++index )
        {
            if ((!fEventInProgress) && (fRunInProgress) && (!fPreEventInProgress))
                {
                    if (ReceivedKassReady()) fPreEventInProgress = true;
                }

            if (fPreEventInProgress)
            {
                PreEventCounter += 1;
//                 printf("preeventcounter is %d\n", PreEventCounter);
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
                        DriveAntenna(PreEventCounter, index, aSignal);

                        PreEventCounter = 0; // reset
                    }
                    tLock.unlock();
                }

            }  // for loop

        // trigger any remaining events in Kassiopeia so that its thread can finish.
        while (fRunInProgress)
        {
            if (fRunInProgress)
            {
            	std::this_thread::sleep_for(std::chrono::milliseconds(100));
            	if (!fEventInProgress)
                  if (ReceivedKassReady())
            	    WakeBeforeEvent();
            }
        }

        Kassiopeia.join();  // wait for Kassiopeia to finish.

        return true;
    }

} /* namespace locust */

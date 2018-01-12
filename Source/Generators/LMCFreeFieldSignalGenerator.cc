/*
 * LMCFreeFieldSignalGenerator.cc
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#include "LMCFreeFieldSignalGenerator.hh"
#include "LMCEventHold.hh"
#include "LMCRunKassiopeia.hh"

#include "logger.hh"
#include <thread>
#include <algorithm>

#include <iostream>
#include <fstream>

#include "LMCGlobalsDeclaration.hh"
#include "LMCHFSSReader.hh"


static std::string gxml_filename = "blank.xml";
const double dx=0.01;
const int NPreEventSamples = 150000;
int fNFDIndex=-1;

namespace locust
{
    LOGGER( lmclog, "FreeFieldSignalGenerator" );

    MT_REGISTER_GENERATOR(FreeFieldSignalGenerator, "freefield-signal");

    FreeFieldSignalGenerator::FreeFieldSignalGenerator( const std::string& aName ) :
            Generator( aName ),
            fWriteNFD(0.),
            fLO_Frequency( 0.)
    {
        fRequiredSignalState = Signal::kTime;
    }

    FreeFieldSignalGenerator::~FreeFieldSignalGenerator()
    {
    }

    bool FreeFieldSignalGenerator::Configure( const scarab::param_node* aParam )
    {
        if( aParam == NULL) return true;
        fLO_Frequency = LO_FREQUENCY;
        if( aParam->has( "lo-frequency" ) )
        {
            fLO_Frequency = aParam->get_value< double >( "lo-frequency" );
        }
        if( aParam->has( "carrier-frequency" ) )
        {
            double fCarrier_Frequency = aParam->get_value< double >( "carrier-frequency" );
            fDecimationFactor = AntiAliasingSetup(fCarrier_Frequency,1.e9);
            LDEBUG( lmclog, "Changing Decimation Factor to: "<< fDecimationFactor);
            fDigitizerTimeStep*= 10. / double(fDecimationFactor);
        }
        if( aParam->has( "xml-filename" ) )
        {
            gxml_filename = aParam->get_value< std::string >( "xml-filename" );
        }

        if( aParam->has( "and-filename" ) )
        {
            fWriteNFD = true;
            fAND_filename = aParam->get_value< std::string >( "and-filename" );
            HFSSReader HFRead;
            HFRead.ParseANDFile(fAND_filename);
            rReceiver=HFRead.GetSurfacePoints();
            NFDFrequencies=HFRead.GetFrequencies();
            fNFD_filename=HFRead.GetNFDFilename();
        }
        else
        {
            //If not, use existing code to generate plane receiver
            HFSSReader HFRead;
            rReceiver = HFRead.GeneratePlane({dx,dx},7);//Argumemts: Size, resolution
            rReceiver = HFRead.RotateShift(rReceiver,{1.,0.,0.},{0.05,0.,0.});//Arguments Normal vector, Position (m)
        }
        PreviousTimes = std::vector<std::pair<int,double> >(rReceiver.size(),{-99.,-99.});

        return true;
    }

    void FreeFieldSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


  static void* KassiopeiaInit()
    {
        //cout << gxml_filename; getchar();
        const std::string & afile = gxml_filename;
    	RunKassiopeia *RunKassiopeia1 = new RunKassiopeia;
    	RunKassiopeia1->Run(afile);
    	delete RunKassiopeia1;
    }



    static void WakeBeforeEvent()
    {
        fPreEventCondition.notify_one();
        return;
    }


    static bool ReceivedKassReady()
    {
        if( !fKassEventReady)
        {
            std::unique_lock< std::mutex >tLock( fKassReadyMutex );
            fKassReadyCondition.wait( tLock );
            printf("LMC Got the fKassReadyCondition signal\n");
        }

        return true;
    }


    //Change decimation factor/ (and therefore the sampling rate) to guarantee no aliasing of signal
    double FreeFieldSignalGenerator::AntiAliasingSetup(double fCarrierFrequency, double fBandwidthFrequency) const
    {
        double fSampleMin = 2. * fBandwidthFrequency;
        int fDecimationRange[2] = {10 , 100};
        int fDecimation = fDecimationRange[0];

        double fSampleFrequency;
        int mFactor;

        while(fDecimation <= fDecimationRange[1])
        {
            fSampleFrequency = fSampleMin *  fDecimation  / fDecimationRange[0];
            mFactor = floor( ( 2. * fCarrierFrequency - fBandwidthFrequency ) / fSampleFrequency);
            if(fSampleFrequency >= ( 2. * fCarrierFrequency + fBandwidthFrequency ) / (mFactor + 1.))
            {
                break;
            }

            ++fDecimation;
        }

        if(fDecimation==fDecimationRange[1])
        {
                LERROR( lmclog, "Cannot find Decimation Factor. Are you sure about your carrier frequency?");
        }

        return fDecimation;
    }

    double m(double angle, double x0, double y0)
    {
    angle = 3.1415926*angle/180.;
    double xprime = cos(angle)*x0 - sin(angle)*y0;
    double yprime = sin(angle)*x0 + cos(angle)*y0;

    double m = yprime/xprime;

    return m;
    }

    double b(double m, double x0, double y0)
    {
    double b = y0 - m*x0;
    return b;
    }

    double directivity(double m1, double m2, double x0, double y0, double x, double y)
    {

    	// cone shaped directivity with gain=1 and half angle OpeningAngle for antenna at x0,y0

    	double directivity = 0.;

    	if (((y < m1*x + b(m1,x0,y0)) && (y > m2*x + b(m2,x0,y0))) |
    	   ((y > m1*x + b(m1,x0,y0)) && (y < m2*x + b(m2,x0,y0))))
    	  {
    	  if (fabs(atan(m1)-atan(m2)) < 1.57)
    	    directivity = 1.;
    	  }
    	else
    	  {
    	  if (fabs(atan(m1)-atan(m2)) > 1.57)
    	    directivity = 1.;
    	  }
    	return directivity;

    }



    void* FreeFieldSignalGenerator::FilterNegativeFrequencies(Signal* aSignal, double *ImaginarySignal) const
    {

/*
        int nwindows = 80;
        int windowsize = 10*aSignal->TimeSize()/nwindows;


        fftw_complex *SignalComplex;
        SignalComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * windowsize );
        fftw_complex *FFTComplex;
        FFTComplex = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * windowsize );

        fftw_plan ForwardPlan;
        ForwardPlan = fftw_plan_dft_1d(windowsize, SignalComplex, FFTComplex, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan ReversePlan;
        ReversePlan = fftw_plan_dft_1d(windowsize, FFTComplex, SignalComplex, FFTW_BACKWARD, FFTW_ESTIMATE);


        for (int nwin = 0; nwin < nwindows; nwin++)
        {
            // Construct complex voltage.

            for( unsigned index = 0; index < windowsize; ++index )
            {
                SignalComplex[index][0] = aLongSignal[ nwin*windowsize + index ];
                SignalComplex[index][1] = ImaginarySignal[ nwin*windowsize + index ];
                //if (index==20000) {printf("signal 20000 is %g\n", aSignal->SignalTime()[index]); getchar();}
            }

            fftw_execute(ForwardPlan);

            //Complex filter to set power at negative frequencies to 0.

            for( unsigned index = windowsize/2; index < windowsize; ++index )
            {
                FFTComplex[index][0] = 0.;
                FFTComplex[index][1] = 0.;
            }

            fftw_execute(ReversePlan);

            double norm = (double)(windowsize);

            for( unsigned index = 0; index < windowsize; ++index )
            {
                // normalize and take the real part of the reverse transform, for digitization.
                //aSignal->SignalTime()[ nwin*windowsize + index ] = SignalComplex[index][0]/norm;
                aLongSignal[ nwin*windowsize + index ] = SignalComplex[index][0]/norm;
                //if (index>=20000) {printf("filtered signal is %g\n", aSignal->SignalTime()[index]); getchar();}
            }

        }


        delete SignalComplex;
        delete FFTComplex;
*/


    }


    void FreeFieldSignalGenerator::NFDWrite() const
    {
            std::ofstream fNFDOutput;
            fNFDOutput.open(fNFD_filename,std::ios::out | std::ios::trunc);

            fNFDOutput << "#Index, X, Y, Z, Ex(real, imag), Ey(real, imag), Ez(real, imag), Hx(real, imag), Hy(real, imag), Hz(real, imag)" <<std::endl;
            
            fNFDOutput << std::fixed;
            fNFDOutput << "Frequencies "<<NFDFrequencies.size()<<std::endl;

            fNFDOutput << std::setprecision(9);

            // Loop over frequencies
            for(int i = 0; i < NFDFrequencies.size() ; ++i)
            {
                fNFDOutput<<std::scientific;
                fNFDOutput <<"Frequency "<< NFDFrequencies[i] <<std::endl;

                for(int j = 0; j < rReceiver.size() ; ++j)
                {

                    fNFDOutput<<std::scientific;
                    fNFDOutput <<j+1<<", "<<rReceiver[j][0]<<", "<<rReceiver[j][1]<<", "<<rReceiver[j][2]<<", "<< NFDElectricFieldFreq[i][j][0][0]<<", "<< NFDElectricFieldFreq[i][j][0][1]<<", "<< NFDElectricFieldFreq[i][j][1][0]<<", "<< NFDElectricFieldFreq[i][j][1][1]<<", "<< NFDElectricFieldFreq[i][j][2][0]<<", "<< NFDElectricFieldFreq[i][j][2][1]<<", "<< NFDMagneticFieldFreq[i][j][0][0]<<", "<< NFDMagneticFieldFreq[i][j][0][1]<<", "<< NFDMagneticFieldFreq[i][j][1][0]<<", "<< NFDMagneticFieldFreq[i][j][1][1]<<", "<< NFDMagneticFieldFreq[i][j][2][0]<<", "<< NFDMagneticFieldFreq[i][j][2][1]<<std::endl;

                }
            }
            fNFDOutput.close();


            //Be 100% sure vector deallocates memory correctly (probably overkill)
            NFDElectricFieldFreq.resize(0);
            NFDMagneticFieldFreq.resize(0);

            NFDElectricFieldFreq.shrink_to_fit();
            NFDMagneticFieldFreq.shrink_to_fit();

    }

    double GetSpaceTimeInterval(const double &aParticleTime, const double &aReceiverTime, const KGeoBag::KThreeVector &aParticlePosition, const KGeoBag::KThreeVector &aReceiverPosition )
    {
        //return pow(aReceiverTime - aParticleTime,2.) - (aReceiverPosition - aParticlePosition).MagnitudeSquared() / pow(KConst::C() , 2.);
        return aReceiverTime - aParticleTime - (aReceiverPosition - aParticlePosition).Magnitude() / KConst::C();
    }

    double GetStepRoot(const locust::Particle aParticle, double aReceiverTime, KGeoBag::KThreeVector aReceiverPosition, double aSpaceTimeInterval, const int aStepOrder = 0)
    {
        double tRetardedTime = aParticle.GetTime(true); //interpolate!!!

        double c=KConst::C();

        if(aStepOrder==0)
        {
            double tCorrection = sqrt(fabs(aSpaceTimeInterval));
            double tSign;
            (aSpaceTimeInterval > 0) ? tSign = 1. : tSign = -1.;
            
            return tRetardedTime + aSpaceTimeInterval;
            //return tRetardedTime + tSign * aSpaceTimeInterval;
        }

        KGeoBag::KThreeVector tNewPosition = aParticle.GetPosition(true);
        KGeoBag::KThreeVector tNewVelocity = aParticle.GetVelocity(true);

        KGeoBag::KThreeVector tReceiverVector = aReceiverPosition - tNewPosition;
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
        KGeoBag::KThreeVector tNewAcceleration = aParticle.GetAcceleration(true);
        double fZeroDoublePrime = 2. * (1. - tNewVelocity.Dot(tNewVelocity)/(c*c)-tNewAcceleration.Dot(tNewPosition-aReceiverPosition)/(c*c));

        if(aStepOrder==2)
        {
            return tRetardedTime-tNewtonRatio* ( 1. + ( tNewtonRatio * fZeroDoublePrime) / ( 2. * fZeroPrime));
        }

        LERROR( lmclog, "Need to put root finding method with order 0-2!" );
        return 0;
    }


    void* FreeFieldSignalGenerator::DriveAntenna(int PreEventCounter, unsigned index, Signal* aSignal) const
    {

        locust::Particle tCurrentParticle = fParticleHistory.back();
        int CurrentIndex;
        HFSSReader HFRead;

        const int nelementsZ = 9;
        double theta = 0.;  // azimuthal angle of each amplifier channel.
        double radius = 0.05; // radius of each amplifier channel in xy plane.
        double OpeningAngle = 15.; //degrees, half angle fake directivity angle.
        double DirectivityFactor = 0.;
        double m1 = 0.;
        double m2 = 0.;
        static double tVoltagePhase[nelementsZ*NCHANNELS] = {0.};
        static double phi_LO = 0.;
        int signalSize = aSignal->TimeSize();

        //Receiver Properties

        double tReceiverTime = t_old;
        KGeoBag::KThreeVector tReceiverPosition;

        double tRetardedTime = 0.; //Retarded time of particle corresponding to when emission occurs, reaching receiver at tReceiverTime
        double tTotalPower=0.;

        double tSpaceTimeInterval=99.;
        double dtRetarded=0;
        double tTolerance=1e-23;
        const double dtStepSize = fabs(fParticleHistory[0].GetTime() - fParticleHistory[1].GetTime());

        const int HistorySize = fParticleHistory.size();

        phi_LO+= 2. * KConst::Pi() * fLO_Frequency * fDigitizerTimeStep;  // this has to happen outside the signal generating loop.

//        printf("temp is %g and phi_LO is %g\n", temp, phi_LO);
//        printf("inferred fLO_Frequency is %g\n", (phi_LO-temp)/2./KConst::Pi()/fDigitizerTimeStep);

       //printf("Size: %d %d\n",HistorySize, fParticleHistory.size());

        //int tAverageIterations=0; //Performance tracker. Count number of iterations to converge....

        for (int ch=0; ch<NCHANNELS; ch++)
        {
        for (int z_position = 4; z_position<5; z_position++) // step through antennas along z
        {
        // position waveguide in space:
        rReceiver = HFRead.GeneratePlane({dx,dx},7);//Argumemts: Size, resolution
        theta = (double)ch*360./NCHANNELS*PI/180.;
        rReceiver = HFRead.RotateShift(rReceiver,{cos(theta),sin(theta),0.},{radius*cos(theta),radius*sin(theta),(double)(z_position-4)*0.01});//Arguments Normal vector, Position (m)
        PreviousTimes = std::vector<std::pair<int,double> >(rReceiver.size(),{-99.,-99.}); // initialize
        tTotalPower = 0.; // initialize
        m1=m(-OpeningAngle*2.,radius*cos(theta), radius*sin(theta));
        m2=m(OpeningAngle*2.,radius*cos(theta), radius*sin(theta));


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

                tRetardedTime = tReceiverTime - (tCurrentParticle.GetPosition() - tReceiverPosition).Magnitude() /  KConst::C();
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

                tRetardedTime = GetStepRoot(tCurrentParticle, tReceiverTime, tReceiverPosition, tSpaceTimeInterval,0);

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


            KGeoBag::KThreeVector tECrossH = tCurrentParticle.CalculateElectricField(rReceiver[i]).Cross(tCurrentParticle.CalculateMagneticField(rReceiver[i]));
            KGeoBag::KThreeVector tDirection = tReceiverPosition - tCurrentParticle.GetPosition(true);

            tTotalPower += dx * dx * tECrossH.Dot(tDirection.Unit()) / rReceiver.size() ;// * (fabs(tCurrentParticle.GetPosition(true).Z())<0.01);

            double tVelZ = tCurrentParticle.GetVelocity(true).Z();
            double tCosTheta =  tVelZ * tDirection.Z() /  tDirection.Magnitude() / fabs(tVelZ);
            double tDopplerFrequency  = tCurrentParticle.GetCyclotronFrequency() / ( 1. - fabs(tVelZ) / KConst::C() * tCosTheta);
            //double tDopplerFrequency  = tCurrentParticle.GetCyclotronFrequency();
            double tTransitTime = (tReceiverPosition - tCurrentParticle.GetPosition(true)).Magnitude() / KConst::C();

            if (tRetardedTime+tTransitTime > fDigitizerTimeStep)  // if the signal has been present for longer than fDigitizerTimeStep
              {
              tVoltagePhase[ch*nelementsZ + z_position]+= tDopplerFrequency * fDigitizerTimeStep / rReceiver.size();
              }
            else  // if this is the first digitizer sample
              {
                tVoltagePhase[ch*nelementsZ + z_position]+=
                		tDopplerFrequency * (tReceiverTime - (tRetardedTime+tTransitTime)) / rReceiver.size();
//                printf("Retarding voltage phase:  fDigitizerTimeStep is %g and tReceiverTime is %g and tRetardedTime is %g\n and t_old is %g\n",
//                		fDigitizerTimeStep, tReceiverTime, tRetardedTime, t_old); getchar();
              }

            /*
            if (i==0)  // check receiver point 0 for each channel.  It should be the same each time.
            {
              printf("rx point 0:  zposition is %d and channel is %d, fcyc is %g and tVelZ is %g, dopplerfreq is %g, costheta is %f\n",
            		  z_position, ch, tDopplerFrequency, tCurrentParticle.GetCyclotronFrequency(), tVelZ, tCosTheta);
              printf("current index is %d\n", CurrentIndex);
            }
            */


            double tMinHFSS=1e-8;
            const int nHFSSBins=2048;


            if( fWriteNFD && (fNFDIndex < nHFSSBins) && (tReceiverTime >= tMinHFSS) )
            {
                KGeoBag::KThreeVector tmpElectricField, tmpMagneticField;
                tmpElectricField = tCurrentParticle.CalculateElectricField(rReceiver[i]);
                tmpMagneticField = tCurrentParticle.CalculateMagneticField(rReceiver[i]);
                
                //If doing the first receiver
                if(!i) ++fNFDIndex;

                for(int j=0;j<NFDFrequencies.size();j++)
                {
                    for(int k=0;k<3;k++)
                    {
                        //Complex factors for Downmixing/ DFT sums
                        double DownConvert[2] = { cos(phi_LO) , - sin(phi_LO) };
                        int kFreq = ( NFDFrequencies[j] - fLO_Frequency ) * nHFSSBins * fDigitizerTimeStep;
                        double DFTFactor[2] = { cos(2. * KConst::Pi() * kFreq * fNFDIndex  / nHFSSBins), -sin(2. * KConst::Pi() * kFreq * fNFDIndex / nHFSSBins ) };

                        NFDElectricFieldFreq[j][i][k][0] += tmpElectricField[k] * ( DownConvert[0] * DFTFactor[0] - DownConvert[1] * DFTFactor[1] ) / nHFSSBins;
                        NFDElectricFieldFreq[j][i][k][1] += tmpElectricField[k] * ( DownConvert[0] * DFTFactor[1] + DownConvert[1] * DFTFactor[0] ) / nHFSSBins;

                        NFDMagneticFieldFreq[j][i][k][0] += tmpMagneticField[k] * ( DownConvert[0] * DFTFactor[0] - DownConvert[1] * DFTFactor[1] ) / nHFSSBins;
                        NFDMagneticFieldFreq[j][i][k][1] += tmpMagneticField[k] * ( DownConvert[0] * DFTFactor[1] + DownConvert[1] * DFTFactor[0] ) / nHFSSBins;
                    }
                }
                //if(fNFDIndex%100==0 && i==0)printf("%d\n",fNFDIndex);

            }
            else if(fWriteNFD && fNFDIndex >=nHFSSBins)
            {
                printf("Done! \n");
            }

        }  // i, rReceiver.size() loop.



        DirectivityFactor = directivity(m1, m2, radius*cos(theta*PI/180.), radius*sin(theta*PI/180.),
        		tCurrentParticle.GetPosition(true).GetX(), tCurrentParticle.GetPosition(true).GetY());

//        if (fabs(tCurrentParticle.GetPosition(true).GetZ()-(double)(z_position-4)*0.01) < 0.005 )
//        {
          aSignal->LongSignalTimeComplex()[ch*signalSize*10 + index][0] += DirectivityFactor * sqrt(50.)*sqrt(tTotalPower) * cos(tVoltagePhase[ch*nelementsZ + z_position] - phi_LO);
          aSignal->LongSignalTimeComplex()[ch*signalSize*10 + index][1] += DirectivityFactor * sqrt(50.)*sqrt(tTotalPower) * sin(tVoltagePhase[ch*nelementsZ + z_position] - phi_LO);
//        }

        } // z_position waveguide element stepping loop.
        } // NCHANNELS loop.

//        printf("tTotalPower at time %g is %g\n\n\n", t_old, tTotalPower); getchar();

        /*
        printf("signal at time %g is %g\n\n\n", t_old, aSignal->LongSignalTimeComplex()[index][0]);
        printf("phi_LO is %g\n", phi_LO);
        printf("tVoltagePhase is %f\n", tVoltagePhase);
        printf("fDigitizerTimeStep is %g\n", fDigitizerTimeStep);getchar();
        */


        //        printf("tTotalPower at time %g is %g\n\n\n", t_old, tTotalPower); getchar();

        t_old += fDigitizerTimeStep;

    }


    //Return iterator of fParticleHistory particle closest to the time we are evaluating
    int FreeFieldSignalGenerator::FindNode(double tNew, double dtStepSize, int kIndexOld) const
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
        //printf("%d %e %e %e\n",tNodeIndex, tNew, fParticleHistory.front().GetTime(), fParticleHistory.back().GetTime());
        //cout<<tNodeIndex<<endl;
        //if(tNodeIndex < 0 || tNodeIndex > (tHistorySize - 1))
        //{
        //    LERROR(lmclog, "OUT OF INDEX SEARCH");
        //}

        return tNodeIndex;
    }

    bool FreeFieldSignalGenerator::DoGenerate( Signal* aSignal ) const
    {

        if(fWriteNFD)
        {
            //Initialize to correct sizes/ all zeros
            NFDElectricFieldFreq.resize(NFDFrequencies.size());
            NFDMagneticFieldFreq.resize(NFDFrequencies.size());
            for(int i=0;i<NFDFrequencies.size();++i)
            {
                NFDElectricFieldFreq[i]=std::vector<std::array<std::array<double,2>, 3> >(rReceiver.size(),{0.});
                NFDMagneticFieldFreq[i]=std::vector<std::array<std::array<double,2>, 3> >(rReceiver.size(),{0.});
            }
        }

        //n samples for event spacing.
        int PreEventCounter = 0;
//        bool fPhaseIISimulation = false;  // this is now a parameter in the xml file.

//        printf("fwritenfd is %d\n", fWriteNFD); getchar();

        std::thread Kassiopeia (KassiopeiaInit);     // spawn new thread
        fRunInProgress = true;

        for( unsigned index = 0; index < 10*aSignal->TimeSize(); ++index )
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

        if(fWriteNFD) NFDWrite();
        
        //FilterNegativeFrequencies(aSignal, ImaginarySignal);
        //delete [] ImaginarySignal;
        
        // trigger any remaining events in Kassiopeia so that its thread can finish.
        while (fRunInProgress)
        {
            if (fRunInProgress)
                if (ReceivedKassReady()) WakeBeforeEvent();
        }

        Kassiopeia.join();  // wait for Kassiopeia to finish.

        return true;
    }

} /* namespace locust */

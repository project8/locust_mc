/*
 * LMCKassSignalGenerator.cc
 *
 *  Created on: Mar 12, 2014
 *      Author: nsoblath
 */

#include "LMCKassSignalGenerator.hh"
#include "LMCEventHold.hh"
#include "LMCRunKassiopeia.hh"
#include "LMCGlobals.hh"

#include "logger.hh"

namespace locust
{
    LOGGER( lmclog, "KassSignalGenerator" );

    MT_REGISTER_GENERATOR(KassSignalGenerator, "kass-signal");

    KassSignalGenerator::KassSignalGenerator( const std::string& aName ) :
            Generator( aName )
    {
        fRequiredSignalState = Signal::kTime;
    }

    KassSignalGenerator::~KassSignalGenerator()
    {
    }

    bool KassSignalGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;

        return true;
    }

    void KassSignalGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    void* KassiopeiaInit(void*)
    {
    	const string & afile = "/home/penny/project8/Project8_by_Devin.xml";
    	RunKassiopeia *RunKassiopeia1 = new RunKassiopeia;
    	RunKassiopeia1->Run(afile);
    	delete RunKassiopeia1;
    }



    bool KassSignalGenerator::DoGenerate( Signal* aSignal ) const
    {

        double dt = 5.e-9; // seconds, this should be coming from RunLengthCalculator or something like it.
        double phi_t = 0.; // antenna voltage phase in radians.
        double phiLO_t = 0.; // voltage phase of LO;
        double fprime = 0.;  // Doppler shifted cyclotron frequency in Hz.

    	pthread_t tid1;
    	pthread_create(&tid1,NULL,KassiopeiaInit ,NULL);


//    	for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        	for( unsigned index = 0; index < 67000; ++index )
        {
            // lock access to the mutex, unless Kassiopeia needs to write to the globals.
                pthread_mutex_lock (&mymutex);
        // stop here and check whether Kassiopeia sent out a tick signaling that the globals need to be digitized.
                pthread_cond_wait(&tick, &mymutex);
//                printf("paused in Locust! zvelocity is %g\n", zvelocity); getchar();
                fprime = fcyc*GammaZ*(1.-zvelocity/2.99792e8);
                phi_t += 2.*PI*fprime*dt;
                phiLO_t += -2.*PI*LO_FREQUENCY*dt;

                aSignal->SignalTime()[ index ] +=
                            pow(LarmorPower,0.5)*(cos(phi_t)*cos(phiLO_t) - sin(phi_t)*sin(phiLO_t));


               printf("Locust says:  index is %d and t is %g and zvelocity is %g and sqrtLarmorPower is %g and fcyc is %.10g and fprime is %g and GammaZ is %f\n",
                           index, t, zvelocity, pow(LarmorPower,0.5), fcyc, fprime, GammaZ);
//               getchar();

                pthread_mutex_unlock (&mymutex);


        }

    	pthread_join(tid1,NULL);  // This makes sure Locust does not just proceed without Kassiopeia.

        return true;
    }

} /* namespace locust */

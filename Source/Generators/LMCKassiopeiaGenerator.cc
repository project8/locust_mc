/*
 * LMCKassiopeiaGenerator.cc
 *
 *  Created on: June 29, 2015
 *      Author: plslocum after nsoblath
 */



#include "LMCKassiopeiaGenerator.hh"
#include "logger.hh"
#include "LMCGlobals.hh"


#include "KSRoot.h"
#include "KMessage.h"
#include "KTextFile.h"




#include "KCommandLineTokenizer.hh"
#include "KXMLTokenizer.hh"
#include "KVariableProcessor.hh"
#include "KIncludeProcessor.hh"
#include "KLoopProcessor.hh"
#include "KConditionProcessor.hh"
#include "KPrintProcessor.hh"
#include "KElementProcessor.hh"
#include "KTagProcessor.hh"
#include "KSSimulation.h"

#ifdef Kommon_USE_ROOT
#include "KFormulaProcessor.hh"
#endif

#include "KSMainMessage.h"
#include "KSToolbox.h"

using namespace Kassiopeia;
using namespace katrin;



using std::string;



namespace locust
{
    LOGGER( lmclog, "KassiopeiaGenerator" );

    MT_REGISTER_GENERATOR(KassiopeiaGenerator, "kassiopeia");

    KassiopeiaGenerator::KassiopeiaGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &KassiopeiaGenerator::DoGenerateTime )
    {
        fRequiredSignalState = Signal::kTime;
    }

    KassiopeiaGenerator::~KassiopeiaGenerator()
    {
    }


// pls hack:  This routine is lifted from Kassiopeia.cxx main().

    void* KassiopeiaInit(void*)
    {

    KCommandLineTokenizer tCommandLine;
//    tCommandLine.ProcessCommandLine( argc, argv );
    tCommandLine.ProcessCommandLine();

    KXMLTokenizer tTokenizer;
    KVariableProcessor tVariableProcessor( tCommandLine.GetVariables() );
    KIncludeProcessor tIncludeProcessor;
    KLoopProcessor tLoopProcessor;
    KConditionProcessor tConditionProcessor;
    KPrintProcessor tPrintProcessor;
    KTagProcessor tTagProcessor;
    KElementProcessor tElementProcessor;

	tVariableProcessor.InsertAfter( &tTokenizer );
	tIncludeProcessor.InsertAfter( &tVariableProcessor );

#ifdef Kommon_USE_ROOT
	KFormulaProcessor tFormulaProcessor;
	tFormulaProcessor.InsertAfter( &tVariableProcessor );
	tIncludeProcessor.InsertAfter( &tFormulaProcessor );
#endif

    tLoopProcessor.InsertAfter( &tIncludeProcessor );
    tConditionProcessor.InsertAfter( &tLoopProcessor );
    tPrintProcessor.InsertAfter( &tConditionProcessor );
    tTagProcessor.InsertAfter( &tPrintProcessor );
    tElementProcessor.InsertAfter( &tTagProcessor );


    mainmsg( eNormal ) << "starting..." << eom;



    KSToolbox::GetInstance();

    KTextFile* tFile;
    for( vector< string >::const_iterator tIter = tCommandLine.GetFiles().begin(); tIter != tCommandLine.GetFiles().end(); tIter++ )
    {
        tFile = new KTextFile();
        tFile->AddToNames( *tIter );
        cout << *tIter << "\n"; // pls hack.  Checking filename.
        tTokenizer.ProcessFile( tFile );  // this line is the program.
        delete tFile;
    }



    KSToolbox::DeleteInstance();

    mainmsg( eNormal ) << "...finished" << eom;

    return 0;

    }


    bool KassiopeiaGenerator::Configure( const scarab::param_node* aParam )
    {

        if( aParam == NULL) return true;
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

        return true;
    }

    void KassiopeiaGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    Signal::State KassiopeiaGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void KassiopeiaGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &KassiopeiaGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &KassiopeiaGenerator::DoGenerateFreq;
        }
        else
        {
            LWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }

    bool KassiopeiaGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }



    void* thread1(void*)
    {
    while(1){
//    printf("Hello!!\n");
    printf("about to lock the mutex\n");
    pthread_mutex_lock (&mymutex);
    testvar = -100000.;
    pthread_mutex_unlock (&mymutex);

    }
    return 0;
    }


    void* thread2(void*)
    {
    while(1){
    printf("How are you?\n");
    pthread_mutex_lock (&mymutex);
    testvar += 1.;
    printf("testvar is %f\n", testvar);
    }
    return 0;
    }


    void* thread3(void*)
    {
    while(1){
    printf("This is the last thread?\n");
    pthread_mutex_lock (&mymutex);
    testvar += 1.;
    printf("testvar is %f\n", testvar);
    }
    return 0;
    }




    int threadtest()
    {
    pthread_t tid1,tid2,tid3;

    pthread_create(&tid1,NULL,thread1 ,NULL);
    pthread_create(&tid2,NULL,thread2,NULL);
    pthread_create(&tid3,NULL,thread3,NULL);
    pthread_join(tid1,NULL);
    pthread_join(tid2,NULL);
    pthread_join(tid3,NULL);
    return 0;
    }


    void* timetrace(void*)
    {
    while(1){
    printf("timetrace says Z is %f\n", Z);
    pthread_mutex_lock (&mymutex);
    }
    return 0;
    }





    bool KassiopeiaGenerator::DoGenerateTime( Signal* aSignal ) const
    {



        Kassiopeia::mainmsg( katrin::eNormal ) << "Hello world from Kassiopeia, inside Locust. " << katrin::eom;  getchar();
        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;
        double dt = RunLengthCalculator1->GetBinWidth(); // seconds
        double phi_t = 0.; // antenna voltage phase in radians.
        double phiLO_t = 0.; // voltage phase of LO;
        double fprime = 0.;  // Doppler shifted cyclotron frequency in Hz.



        pthread_t tid1;
        pthread_create(&tid1,NULL,KassiopeiaInit ,NULL);


//printf("locust is going to wait right here.\n"); getchar();

//        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )

            for( unsigned index = 0; index < 45; ++index )
        {

        	// lock access to the mutex, unless Kassiopeia needs to write to the globals.
            pthread_mutex_lock (&mymutex);
            // stop here and check whether Kassiopeia sent out a tick signaling that the globals need to be digitized.
            pthread_cond_wait(&tick, &mymutex);
            printf("Locust says:  index is %d and t is %g and sqrtLarmorPower is %g\n", index, t, pow(LarmorPower, 0.5));

            fprime = fcyc*GammaZ*(1.-zvelocity/2.99792e8);
            phi_t += 2.*PI*fprime*dt;
            phiLO_t += -2.*PI*LO_FREQUENCY*dt;

            aSignal->SignalTime()[ index ] +=
            		pow(LarmorPower,0.5)*(cos(phi_t)*cos(phiLO_t) - sin(phi_t)*sin(phiLO_t));

            pthread_mutex_unlock (&mymutex);

 //           if (index<10) printf("signal %d is %g\n", index, aSignal->SignalTime(index));



        }



        pthread_join(tid1,NULL);  // This makes sure Locust does not just proceed without Kassiopeia.
        delete RunLengthCalculator1;


        return true;
    }

    bool KassiopeiaGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
         return true;
    }

} /* namespace locust */

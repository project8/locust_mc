/*
 * LMCKassiopeiaGenerator.cc
 *
 *  Created on: June 29, 2015
 *      Author: plslocum after nsoblath
 */




#include "LMCKassiopeiaGenerator.hh"
#include "../Core/LMCLogger.hh"




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
    LMCLOGGER( lmclog, "KassiopeiaGenerator" );

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
    void KassiopeiaGenerator::KassiopeiaInit()
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
        cout << *tIter << "\n"; // pls hack.  Check filename.
        tTokenizer.ProcessFile( tFile );
        delete tFile;
    }


    KSToolbox::DeleteInstance();

    mainmsg( eNormal ) << "...finished" << eom;




    return;


    }


    bool KassiopeiaGenerator::Configure( const ParamNode* aParam )
    {

        KassiopeiaInit();

        if( aParam == NULL) return true;
        if( aParam->Has( "domain" ) )
        {
            string domain = aParam->GetValue( "domain" );
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                LMCDEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                LMCERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
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
            LMCWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }

    bool KassiopeiaGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool KassiopeiaGenerator::DoGenerateTime( Signal* aSignal ) const
    {


//        Kassiopeia::mainmsg( katrin::eNormal ) << "Hello world from Kassiopeia, inside Locust. " << katrin::eom;  getchar();
        return true;
    }

    bool KassiopeiaGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
         return true;
    }

} /* namespace locust */

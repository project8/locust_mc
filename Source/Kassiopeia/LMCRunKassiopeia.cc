/*
 * LMCRunKassiopeia.cc
 *
 *  Created on: Mar 10, 2016
 *      Author: nsoblath
 */

#include "LMCRunKassiopeia.hh"


#include <KSSimulation.h>
#include <KSRoot.h>
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

//#ifdef Kassiopeia_USE_ROOT
#include "KFormulaProcessor.hh"
//#endif

#include "KSMainMessage.h"
#include "KToolbox.h"




using namespace Kassiopeia;
using namespace katrin;
//using namespace std;


namespace locust
{

    RunKassiopeia::RunKassiopeia()/*:
    fTokenizer( new KXMLTokenizer() ),
    fCommandLineTokenizer( new KCommandLineTokenizer() ),
    fVariableProcessor( new KVariableProcessor(fCommandLineTokenizer->GetVariables()) ),
    fIncludeProcessor( new KIncludeProcessor() ),
    fFormulaProcessor( new KFormulaProcessor() ),
    fLoopProcessor( new KLoopProcessor() ),
    fConditionProcessor( new KConditionProcessor() ),
    fPrintProcessor( new KPrintProcessor() ),
    fTagProcessor( new KTagProcessor() ),
    fElementProcessor( new KElementProcessor() )*/
    {
/*
    	char* dummy_args[] = { "dummyname", "/home/penny/locust_mc/Config/Project8Phase2_WithRoot.xml", NULL};
        fCommandLineTokenizer->ProcessCommandLine(2, dummy_args);
        fCommandLineTokenizer->GetVariables();




    	fVariableProcessor->InsertAfter( fTokenizer );
        fIncludeProcessor->InsertAfter( fVariableProcessor );


//    #ifdef Kommon_USE_ROOT
        fFormulaProcessor->InsertAfter( fVariableProcessor );
        fIncludeProcessor->InsertAfter( fFormulaProcessor );
//    #endif



        fLoopProcessor->InsertAfter( fIncludeProcessor );
        fConditionProcessor->InsertAfter( fLoopProcessor );
        fPrintProcessor->InsertAfter( fConditionProcessor );
        fTagProcessor->InsertAfter( fPrintProcessor );
        fElementProcessor->InsertAfter( fTagProcessor );
*/
    }


    RunKassiopeia::~RunKassiopeia()
    {
    	/*
        delete fTokenizer;
        delete fVariableProcessor;
        delete fIncludeProcessor;
        delete fLoopProcessor;
        delete fConditionProcessor;
        delete fPrintProcessor;
        delete fTagProcessor;
        delete fElementProcessor;

//#ifdef Kommon_USE_ROOT
        delete fFormulaProcessor;
//#endif
  */

        KToolbox::GetInstance().Clear();
    }

<<<<<<< HEAD
    void RunKassiopeia::SetVariableMap( const std::map< std::string, std::string >& aMap )
    {
        fVariableProcessor->SetExternalMap( aMap );
        return;
    }
=======

>>>>>>> my-temporary-work

    int RunKassiopeia::Run( const std::vector< std::string >& aFiles )
    {

<<<<<<< HEAD
        KToolbox::GetInstance();

        KTextFile* tFile;
        for( std::vector< std::string >::const_iterator tIter = aFiles.begin(); tIter != aFiles.end(); tIter++ )
=======

//    cout << "file vector is \n";
//    cout << aFiles[0]; getchar();

    	char* dummy_args[] = { "dummyname", "/home/penny/Kassiopeia_3.3.2/Kassiopeia/XML/Examples/QuadrupoleTrapSimulation.xml", NULL};

        KCommandLineTokenizer tCommandLine;
        tCommandLine.ProcessCommandLine( 2, dummy_args );

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

//#ifdef Kassiopeia_USE_ROOT
	KFormulaProcessor tFormulaProcessor;
	tFormulaProcessor.InsertAfter( &tVariableProcessor );
	tIncludeProcessor.InsertAfter( &tFormulaProcessor );
//#endif


    tLoopProcessor.InsertAfter( &tIncludeProcessor );
    tConditionProcessor.InsertAfter( &tLoopProcessor );
    tPrintProcessor.InsertAfter( &tConditionProcessor );
    tTagProcessor.InsertAfter( &tPrintProcessor );

    tElementProcessor.InsertAfter( &tTagProcessor );



       mainmsg( eNormal ) << "starting ..." << eom;
       KToolbox::GetInstance();



        KTextFile* tFile;



        for( std::vector< std::string >::const_iterator tIter = aFiles.begin(); tIter != aFiles.end(); tIter++ )

>>>>>>> my-temporary-work
        {
            tFile = new KTextFile();
            tFile->AddToNames( *tIter );
            tTokenizer.ProcessFile( tFile );
            delete tFile;
        }

        mainmsg( eNormal ) << "... finished" << eom;


        return 0;
    }

    int RunKassiopeia::Run( const std::string& aFile )
    {
<<<<<<< HEAD
=======
    	
>>>>>>> my-temporary-work
        std::vector< std::string > tFileVec( 1 );
        tFileVec[ 0 ] = aFile;
        return Run( tFileVec );
        
    }

} /* namespace locust */

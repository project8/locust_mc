/*
 * LMCRunKassiopeia.cc
 *
 *  Created on: Mar 10, 2016
 *      Author: nsoblath
 */

#include "LMCRunKassiopeia.hh"

#include "KSSimulation.h"
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

    RunKassiopeia::RunKassiopeia()
    {

    }


    RunKassiopeia::~RunKassiopeia()
    {
        KToolbox::GetInstance().Clear();
    }


    int RunKassiopeia::Run( const std::vector< std::string >& aFiles, kl_interface_ptr_t aInterface )
    {


        //    cout << "file vector is \n";
        //    cout << aFiles[0].c_str();

        char* dummy_args[] = { (char*)"dummyname", (char*)aFiles[0].c_str(), NULL};

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

        KTextFile* tFile;

        mainmsg( eNormal ) << "starting ..." << eom;

        for( std::vector< std::string >::const_iterator tIter = aFiles.begin(); tIter != aFiles.end(); tIter++ )
        {
            tFile = new KTextFile();
            tFile->AddToNames( *tIter );
            tTokenizer.ProcessFile( tFile );
            delete tFile;
        }

        mainmsg( eNormal ) << "... finished" << eom;

        //	tTokenizer.~KXMLTokenizer();

        KToolbox::GetInstance().Clear();

        return 0;
    }

    int RunKassiopeia::Run( const std::string& aFile, kl_interface_ptr_t aInterface )
    {
        std::vector< std::string > tFileVec( 1 );
        tFileVec[ 0 ] = aFile;
        return Run( tFileVec, aInterface );

    }

} /* namespace locust */


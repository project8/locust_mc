/*
 * LMCRunKassiopeia.cc
 *
 *  Created on: Mar 10, 2016
 *      Author: nsoblath
 */

#include "LMCRunKassiopeia.hh"

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
#include "KSaveSettingsProcessor.hh"
#endif

#include "KSMainMessage.h"
#include "KSRoot.h"
#include "KSToolbox.h"

using namespace katrin;
using namespace Kassiopeia;

namespace locust
{

    RunKassiopeia::RunKassiopeia( KCommandLineTokenizer& aCommandLine )
    {
        KXMLTokenizer tTokenizer;
        KVariableProcessor tVariableProcessor( aCommandLine.GetVariables() );
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

    #ifdef Kommon_USE_ROOT
        KSaveSettingsProcessor tSSProcessor;
        tSSProcessor.InsertAfter( &tPrintProcessor );
        tTagProcessor.InsertAfter( &tSSProcessor );
    #endif

        tElementProcessor.InsertAfter( &tTagProcessor );

        mainmsg( eNormal ) << "starting..." << eom;

        KSToolbox::GetInstance();

        KTextFile* tFile;
        for( vector< string >::const_iterator tIter = aCommandLine.GetFiles().begin(); tIter != aCommandLine.GetFiles().end(); tIter++ )
        {
            tFile = new KTextFile();
            tFile->AddToNames( *tIter );
            tTokenizer.ProcessFile( tFile );
            delete tFile;
        }
    }

    RunKassiopeia::~RunKassiopeia()
    {
        KSToolbox::DeleteInstance();
    }

    int RunKassiopeia::Run()
    {
        KSToolbox* tToolbox = KSToolbox::GetInstance();

        KSRoot* tRoot = tToolbox->GetObjectAs< KSRoot >( "KSRoot" );

        if( tRoot == NULL )
        {
            mainmsg( eError ) << "Did not get a KSRoot object" << eom;
            return -1;
        }

        return 0;
    }

} /* namespace locust */

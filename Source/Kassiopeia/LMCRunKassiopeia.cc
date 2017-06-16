/*
 * LMCRunKassiopeia.cc
 *
 *  Created on: Mar 10, 2016
 *      Author: nsoblath
 */

#include "LMCRunKassiopeia.hh"

#include "KMessage.h"
#include "KTextFile.h"

#include "KXMLTokenizer.hh"
#include "KVariableProcessor.hh"
#include "KIncludeProcessor.hh"
#include "KLoopProcessor.hh"
#include "KConditionProcessor.hh"
#include "KPrintProcessor.hh"
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

    RunKassiopeia::RunKassiopeia() :
            fTokenizer( new KXMLTokenizer() ),
            fVariableProcessor( new KVariableProcessor() ),
            fIncludeProcessor( new KIncludeProcessor() ),
            fLoopProcessor( new KLoopProcessor() ),
            fConditionProcessor( new KConditionProcessor() ),
            fPrintProcessor( new KPrintProcessor() ),
            fTagProcessor( new KTagProcessor() ),
            fElementProcessor( new KElementProcessor() )
    {
#ifdef Kommon_USE_ROOT
        fFormulaProcessor = new KFormulaProcessor();
        fSSProcessor = new KSaveSettingsProcessor();
#endif

        fVariableProcessor->InsertAfter( fTokenizer );
        fIncludeProcessor->InsertAfter( fVariableProcessor );

    #ifdef Kommon_USE_ROOT
        fFormulaProcessor->InsertAfter( fVariableProcessor );
        fIncludeProcessor->InsertAfter( fFormulaProcessor );
    #endif

        fLoopProcessor->InsertAfter( fIncludeProcessor );
        fConditionProcessor->InsertAfter( fLoopProcessor );
        fPrintProcessor->InsertAfter( fConditionProcessor );
        fTagProcessor->InsertAfter( fPrintProcessor );

    #ifdef Kommon_USE_ROOT
        fSSProcessor->InsertAfter( fPrintProcessor );
        fTagProcessor->InsertAfter( fSSProcessor );
    #endif

        fElementProcessor->InsertAfter( fTagProcessor );
    }

    RunKassiopeia::~RunKassiopeia()
    {
        delete fTokenizer;
        delete fVariableProcessor;
        delete fIncludeProcessor;
        delete fLoopProcessor;
        delete fConditionProcessor;
        delete fPrintProcessor;
        delete fTagProcessor;
        delete fElementProcessor;

#ifdef Kommon_USE_ROOT
        delete fFormulaProcessor;
        delete fSSProcessor;
#endif

        KSToolbox::DeleteInstance();
    }

    void RunKassiopeia::SetVariableMap( const map< string, string >& aMap )
    {
        fVariableProcessor->SetExternalMap( aMap );
        return;
    }

    int RunKassiopeia::Run( const std::vector< std::string >& aFiles )
    {
        mainmsg( eNormal ) << "starting..." << eom;

        KSToolbox::GetInstance();

        KTextFile* tFile;
        for( vector< string >::const_iterator tIter = aFiles.begin(); tIter != aFiles.end(); tIter++ )
        {
            tFile = new KTextFile();
            tFile->AddToNames( *tIter );
            fTokenizer->ProcessFile( tFile );
            delete tFile;
        }

        return 0;
    }

    int RunKassiopeia::Run( const std::string& aFile )
    {
        vector< string > tFileVec( 1 );
        tFileVec[ 0 ] = aFile;
        return Run( tFileVec );
    }

} /* namespace locust */

/*
 * LMCEggWriter.cc
 *
 *  Created on: Feb 6, 2014
 *      Author: nsoblath
 */

#include "LMCEggWriter.hh"

#include "LMCLogger.hh"
#include "LMCParam.hh"

namespace locust
{
    LMCLOGGER( lmclog, "EggWriter" );

    EggWriter::EggWriter() :
            fState( kClosed ),
            fFilename( "locust_mc.egg" ),
            fDate(),
            fDescription(),
            fRunType( monarch::sRunTypeOther ),
            fBitDepth( 8 ),
            fDataTypeSize( 1 ),
            fAcquisitionRate( 1. ),
            fDuration( 1. ),
            fRecordSize( 1 ),
            fMonarch( NULL ),
            fRecord( NULL )
    {
    }

    EggWriter::~EggWriter()
    {
    }

    bool EggWriter::Configure( const ParamNode* aNode )
    {
        if( fState != kClosed )
        {
            LMCERROR( lmclog, "Cannot configure the writer while a file is open" );
            return false;
        }

        if( aNode == NULL ) return true;

        SetFilename( aNode->get_value( "egg-filename", fFilename ) );
        SetDate( aNode->get_value( "date", fDate ) );
        SetDescription( aNode->get_value( "description", fDescription ) );
        SetRunType( aNode->get_value< monarch::RunType >( "run-type", fRunType ) );
        if( ! SetBitDepth( aNode->get_value< unsigned >( "bit-depth", fBitDepth ) ) )
        {
            LMCERROR( lmclog, "Error setting the bit depth" );
            return false;
        }

        return true;
    }

    bool EggWriter::PrepareEgg()
    {
        if( fState != kClosed )
        {
            LMCERROR( lmclog, "Egg preparation cannot begin while a file is open" );
            return false;
        }

        fMonarch = monarch::Monarch::OpenForWriting( fFilename );

        monarch::MonarchHeader* header = fMonarch->GetHeader();

        // configurable items
        header->SetFilename( fFilename );
        header->SetTimestamp( fDate );
        header->SetDescription( fDescription );
        header->SetRunType( fRunType );
        header->SetDataTypeSize( fDataTypeSize );
        header->SetAcquisitionRate( fAcquisitionRate );
        header->SetRunDuration( fDuration );
        header->SetRecordSize( fRecordSize );

        // non-configurable items
        header->SetAcquisitionMode( monarch::sOneChannel );
        header->SetRunSource( monarch::sSourceSimulation );
        header->SetFormatMode( monarch::sFormatSingle );

        fRecord = fMonarch->GetRecordSeparateOne();

        fState = kPrepared;
        return true;
    }

    bool EggWriter::WriteRecord( const Signal* aSignal )
    {
        if( fState != kWriting || fState != kPrepared )
        {
            LMCERROR( lmclog, "Egg file must be opened before writing records" );
            return false;
        }
        fState = kWriting;



        return fMonarch->WriteRecord();
    }

    bool EggWriter::FinalizeEgg()
    {
        if( fState != kWriting || fState != kPrepared )
        {
            LMCERROR( lmclog, "Egg file must be opened before finalizing" );
            return false;
        }

        fMonarch->Close();
        delete fMonarch;

        fState = kClosed;

        return true;
    }

    monarch::Monarch* EggWriter::GetMonarch() const
    {
        return fMonarch;
    }

    const std::string& EggWriter::GetFilename() const
    {
        return fFilename;
    }

    void EggWriter::SetFilename( const std::string& filename )
    {
        fFilename = filename;
        return;
    }

    const std::string& EggWriter::GetDate() const
    {
        return fDate;
    }

    void EggWriter::SetDate( const std::string& date )
    {
        fDate = date;
        return;
    }

    const std::string& EggWriter::GetDescription() const
    {
        return fDescription;
    }

    void EggWriter::SetDescription( const std::string& desc )
    {
        fDescription = desc;
        return;
    }

    monarch::RunType EggWriter::GetRunType() const
    {
        return fRunType;
    }

    void EggWriter::SetRunType( monarch::RunType runType )
    {
        fRunType = runType;
        return;
    }

    unsigned EggWriter::GetBitDepth() const
    {
        return fBitDepth;
    }

    bool EggWriter::SetBitDepth( unsigned bitDepth )
    {
        if( fBitDepth <= 8 )
        {
            fDataTypeSize = 1;
        }
        else if( fBitDepth <= 16 )
        {
            fDataTypeSize = 2;
        }
        else if( fBitDepth <= 32 )
        {
            fDataTypeSize = 4;
        }
        else if( fBitDepth <= 64 )
        {
            fDataTypeSize = 8;
        }
        else
        {
            LMCERROR( lmclog, "Bit depth must be 64 bits or less" );
            return false;
        }
        fBitDepth = bitDepth;
        return true;
    }

    unsigned EggWriter::GetDataTypeSize() const
    {
        return fDataTypeSize;
    }

    double EggWriter::GetAcquisitionRate() const
    {
        return fAcquisitionRate;
    }

    void EggWriter::SetAcquisitionRate( double rate )
    {
        fAcquisitionRate = rate;
        return;
    }

    double EggWriter::GetDuration() const
    {
        return fDuration;
    }

    void EggWriter::SetDuration( double duration )
    {
        fDuration = duration;
        return;
    }

    unsigned EggWriter::GetRecordSize() const
    {
        return fRecordSize;
    }

    void EggWriter::SetRecordSize( unsigned size )
    {
        fRecordSize = size;
        return;
    }


} /* namespace locust */

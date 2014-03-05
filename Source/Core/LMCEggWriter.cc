/*
 * LMCEggWriter.cc
 *
 *  Created on: Feb 6, 2014
 *      Author: nsoblath
 */

#include "LMCEggWriter.hh"

#include "LMCLogger.hh"
#include "LMCParam.hh"
#include "LMCSignal.hh"

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
            fVoltageMin( -0.25 ),
            fVoltageRange( 0.5 ),
            fAcquisitionRate( 1. ),
            fDuration( 1. ),
            fRecordSize( 1 ),
            fRecordLength( 1 ),
            fMonarch( NULL ),
            fRecord( NULL ),
            fAcquisitionId( 0 ),
            fRecordCounter( 0 ),
            fRecordTime( 0. ),
            fRecordNBytes( 0 )
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

        SetFilename( aNode->GetValue( "egg-filename", fFilename ) );
        SetDate( aNode->GetValue( "date", fDate ) );
        SetDescription( aNode->GetValue( "description", fDescription ) );
        SetRunType( aNode->GetValue< monarch::RunType >( "run-type", fRunType ) );

        return true;
    }

    bool EggWriter::PrepareEgg()
    {
        if( fState != kClosed )
        {
            LMCERROR( lmclog, "Egg preparation cannot begin while a file is open" );
            return false;
        }

        LMCDEBUG( lmclog, "Preparing egg file <" << fFilename << ">" );

        fMonarch = monarch::Monarch::OpenForWriting( fFilename );

        monarch::MonarchHeader* header = fMonarch->GetHeader();

        // configurable items
        header->SetFilename( fFilename );
        if( fDate.empty() )
        {
            char timestamp [256];
            get_time_absolute_str( timestamp );
            header->SetTimestamp( string( timestamp ) );
        }
        else
        {
            header->SetTimestamp( fDate );
        }
        header->SetDescription( fDescription );
        header->SetRunType( fRunType );
        header->SetAcquisitionRate( fAcquisitionRate );
        header->SetRunDuration( fDuration );
        header->SetRecordSize( fRecordSize );

        // non-configurable items
        header->SetAcquisitionMode( monarch::sOneChannel );
        header->SetRunSource( monarch::sSourceSimulation );
        header->SetFormatMode( monarch::sFormatSingle );

        header->SetBitDepth( fBitDepth );
        header->SetDataTypeSize( fDataTypeSize );
        header->SetVoltageMin( fVoltageMin );
        header->SetVoltageRange( fVoltageRange );

        fMonarch->WriteHeader();

        fMonarch->SetInterface( monarch::sInterfaceSeparate );
        fRecord = fMonarch->GetRecordSeparateOne();

        fAcquisitionId = 0;
        fRecordCounter = 0;
        fRecordTime = 0;
        fRecordLength = ( double )fRecordSize / ( 1.e3 * fAcquisitionRate );

        fRecordNBytes = header->GetDataTypeSize() * fRecordSize;

        fState = kPrepared;

        LMCINFO( lmclog, "Egg file <" << fFilename << "> is ready for records" );

        return true;
    }

    bool EggWriter::WriteRecord( const Signal* aSignal )
    {
        if( fState != kWriting && fState != kPrepared )
        {
            LMCERROR( lmclog, "Egg file must be opened before writing records" );
            return false;
        }
        fState = kWriting;

        LMCDEBUG( lmclog, "Writing record " << fRecordCounter    );

        fRecord->fAcquisitionId = fAcquisitionId;
        fRecord->fRecordId = fRecordCounter;
        fRecord->fTime = fRecordTime;

        if( aSignal->GetState() != Signal::kDigital )
        {
            LMCERROR( lmclog, "Signal is not digitized (state = " << aSignal->GetState() << "); no record was written" );
            return false;
        }
        ::memcpy( fRecord->fData, reinterpret_cast< const byte_type* >( aSignal->SignalDigital() ), fRecordNBytes );

        ++fRecordCounter;
        fRecordTime += fRecordLength;

        return fMonarch->WriteRecord();
    }

    bool EggWriter::FinalizeEgg()
    {
        if( fState != kWriting && fState != kPrepared )
        {
            LMCERROR( lmclog, "Egg file must be opened before finalizing" );
            return false;
        }

        LMCDEBUG( lmclog, "Closing egg file" );

        fMonarch->Close();
        delete fMonarch;

        fState = kClosed;

        LMCINFO( lmclog, "Egg file closed" );

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

    double EggWriter::GetVoltageMin() const
    {
        return fVoltageMin;
    }

    void EggWriter::SetVoltageMin( double vMin )
    {
        fVoltageMin = vMin;
        return;
    }

    double EggWriter::GetVoltageRange() const
    {
        return fVoltageRange;
    }

    void EggWriter::SetVoltageRange( double vRange )
    {
        fVoltageRange = vRange;
        return;
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

    void EggWriter::IncrementAcquisitionId()
    {
        ++fAcquisitionId;
        return;
    }


} /* namespace locust */

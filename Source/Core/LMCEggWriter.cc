/*
 * LMCEggWriter.cc
 *
 *  Created on: Feb 6, 2014
 *      Author: nsoblath
 */

#include "LMCEggWriter.hh"

#include "LMCDigitizer.hh"
#include "logger.hh"
#include "param.hh"
#include "LMCRunLengthCalculator.hh"
#include "LMCSignal.hh"
#include "LMCSimulationController.hh"

#include "time.hh"

#include <vector>

namespace locust
{
    LOGGER( lmclog, "EggWriter" );

    EggWriter::EggWriter() :
            f_filename( "/home/hep/heeger/ps48/locust_mc/cbuild/bin/locust_mc_1.egg" ),
            f_date(),
            f_description(),
            f_record_length( 0 ),
            f_record_id( 0 ),
            f_record_time( 0 ),
            f_record_n_bytes( 0 ),
            f_monarch( NULL ),
            f_stream( NULL ),
            f_record( NULL ),
            f_state( kClosed )
    {
    }

    EggWriter::~EggWriter()
    {
        if( f_monarch != NULL )
        {
            FinalizeEgg();
        }
    }

    bool EggWriter::Configure( const scarab::param_node* aNode )
    {
        if( f_state != kClosed )
        {
            LERROR( lmclog, "Cannot configure the writer while a file is open" );
            return false;
        }

        if( aNode == NULL ) return true;

        f_filename = aNode->get_value( "egg-filename", f_filename );
        f_date = aNode->get_value( "date", f_date );
        f_description = aNode->get_value( "description", f_description );

        return true;
    }

    bool EggWriter::PrepareEgg( const RunLengthCalculator* a_rlc, const Digitizer* a_digitizer )
    {

    	bool IQStream = true; // move this into json file as a parameter.

    	SimulationController SimulationController1;
        const unsigned nchannels = SimulationController1.GetNChannels();


        if( f_state != kClosed )
        {
            LERROR( lmclog, "Egg preparation cannot begin while a file is open" );
            return false;
        }

        LDEBUG( lmclog, "Preparing egg file <" << f_filename << ">" );

        f_monarch = monarch3::Monarch3::OpenForWriting( f_filename );

        monarch3::M3Header* header = f_monarch->GetHeader();

        // configurable items
        header->SetFilename( f_filename );
        if( f_date.empty() )
        {
            header->SetTimestamp( scarab::get_absolute_time_string()  );
        }
        else
        {
            header->SetTimestamp( f_date );
        }
        header->SetDescription( f_description );
        header->SetRunDuration( a_rlc->GetDuration() );

        unsigned t_data_type_size = 1;
        bool t_signed_vals = false;
        unsigned t_bit_depth = 8;
        bool t_bits_right_aligned = false;
        if( a_digitizer != NULL )
        {
            t_data_type_size = a_digitizer->DigitizerParams().data_type_size;
            t_signed_vals = a_digitizer->GetADCValuesSigned();
            t_bit_depth = a_digitizer->DigitizerParams().bit_depth;
            t_bits_right_aligned = a_digitizer->DigitizerParams().bits_right_aligned;
        }

        std::vector< unsigned > t_chan_vec;
        uint32_t t_stream_id;

        if (!IQStream)
        {
            t_stream_id = header->AddStream( "locust_mc",
                a_rlc->GetAcquisitionRate(), a_rlc->GetRecordSize(), 1,
                t_data_type_size, t_signed_vals,
                t_bit_depth, t_bits_right_aligned,
                &t_chan_vec );
        }
        else
        {
            t_stream_id = header->AddStream( "locust_mc", nchannels, 1,
                a_rlc->GetAcquisitionRate(), a_rlc->GetRecordSize(), a_rlc->GetSampleSize(),
                t_data_type_size, t_signed_vals,
                t_bit_depth, t_bits_right_aligned,
                &t_chan_vec );
        }


        for( std::vector< unsigned >::const_iterator it = t_chan_vec.begin(); it != t_chan_vec.end(); ++it )
        {
            header->GetChannelHeaders()[ *it ].SetVoltageOffset( a_digitizer->DigitizerParams().v_offset );
            header->GetChannelHeaders()[ *it ].SetVoltageRange( a_digitizer->DigitizerParams().v_range );
//            header->GetChannelHeaders()[ *it ].SetDACGain( 1. );
            header->GetChannelHeaders()[ *it ].SetDACGain( a_digitizer->DigitizerParams().dac_gain );
        }

        f_monarch->WriteHeader();

        f_stream = f_monarch->GetStream( t_stream_id );
        f_record = f_stream->GetStreamRecord();

        f_record_id = 0;
        f_record_time = 0;
        f_record_length = ( double )a_rlc->GetRecordSize() / ( 1.e-3 * a_rlc->GetAcquisitionRate() ); // in ns
        f_record_n_bytes = nchannels * a_digitizer->DigitizerParams().data_type_size * a_rlc->GetRecordSize() * a_rlc->GetSampleSize();


        f_state = kPrepared;

        LINFO( lmclog, "Egg file <" << f_filename << "> is ready for records" );

        return true;
    }

    bool EggWriter::WriteRecord( const Signal* aSignal )
    {


        static bool t_is_new_acq = true;


        if( f_state != kWriting && f_state != kPrepared )
        {
            LERROR( lmclog, "Egg file must be opened before writing records" );
            return false;
        }


        f_state = kWriting;

        LDEBUG( lmclog, "Writing record " << f_record_id );



        f_record->SetRecordId( f_record_id );

        f_record->SetTime( f_record_time );



        if( aSignal->GetState() != Signal::kDigital )
        {
            LERROR( lmclog, "Signal is not digitized (state = " << aSignal->GetState() << "); no record was written" );
            return false;
        }



        if( aSignal->GetDigitalIsSigned() ) ::memcpy( f_record->GetData(), reinterpret_cast< const monarch3::byte_type* >( aSignal->SignalDigitalS() ), f_record_n_bytes );
        else ::memcpy( f_record->GetData(), reinterpret_cast< const monarch3::byte_type* >( aSignal->SignalDigitalUS() ), f_record_n_bytes );

        ++f_record_id;
        f_record_time += f_record_length;


        bool t_return = f_stream->WriteRecord( true ); // pls had to edit false to true
        t_is_new_acq = false;
        return t_return;
    }

    bool EggWriter::FinalizeEgg()
    {
        if( f_state != kWriting && f_state != kPrepared )
        {
            LERROR( lmclog, "Egg file must be opened before finalizing" );
            return false;
        }

        if( f_stream != NULL )
        {
            f_stream->Close();
            f_record = NULL;
            f_stream = NULL;
        }

        LDEBUG( lmclog, "Closing egg file" );

        f_monarch->FinishWriting();
        delete f_monarch;
        f_monarch = NULL;

        f_state = kClosed;

        LINFO( lmclog, "Egg file closed" );

        return true;
    }

} /* namespace locust */

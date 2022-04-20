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

#include "time.hh"

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <stdlib.h>

#include <vector>

namespace locust
{
    LOGGER( lmclog, "EggWriter" );

    EggWriter::EggWriter() :
            f_filename( "locust_mc.egg" ),
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

    bool EggWriter::Configure( const scarab::param_node& aNode )
    {
        if( f_state != kClosed )
        {
            LERROR( lmclog, "Cannot configure the writer while a file is open" );
            return false;
        }

        f_filename = aNode.get_value( "egg-filename", f_filename );
        f_date = aNode.get_value( "date", f_date );
        f_description = aNode.get_value( "description", f_description );

/*        if (f_filename=="locust_mc.egg") // placeholder egg file.
          {
	  const char *homedir;
          char fname[200];
	  if ((homedir = getenv("HOME")) == NULL) 
	    {
	    homedir = getpwuid(getuid())->pw_dir;
            }	              
          int n = sprintf(fname, "locust_mc.egg", homedir);
          f_filename = fname;
          }
*/
        return true;
    }

    bool EggWriter::PrepareEgg( const RunLengthCalculator* a_rlc, const Digitizer* a_digitizer )
    {


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
            header->SetTimestamp( scarab::get_formatted_now()  );
        }
        else
        {
            header->SetTimestamp( f_date );
        }
        header->SetDescription( f_description );
        header->SetRunDuration( a_rlc->GetDuration() );

        unsigned t_data_type_size = 1;
        //bool t_signed_vals = false;
        uint32_t t_data_format = 0;
        unsigned t_bit_depth = 8;
        bool t_bits_right_aligned = false;
        unsigned n_channels = a_rlc->GetNChannels();
        double t_v_range = 1.0;
        double t_v_offset = 0.0;
        double t_dac_gain = t_v_range/(1<<t_bit_depth);
        if( a_digitizer )
        {
            t_data_type_size = a_digitizer->DigitizerParams().data_type_size;
            t_data_format = a_digitizer->GetADCValuesSigned();
            t_bit_depth = a_digitizer->DigitizerParams().bit_depth;
            t_bits_right_aligned = a_digitizer->DigitizerParams().bits_right_aligned;
            t_v_range = a_digitizer->DigitizerParams().v_range;
            t_v_offset = a_digitizer->DigitizerParams().v_offset;
            t_dac_gain = a_digitizer->DigitizerParams().dac_gain;
        } else {
			LDEBUG( lmclog, "Running without digitizer, data saved in float format");
			t_data_format = 2;
			t_data_type_size = sizeof(double);
		}

        std::vector< unsigned > t_chan_vec;
        uint32_t t_stream_id;

            t_stream_id = header->AddStream( "locust_mc", n_channels, 1,
                a_rlc->GetAcquisitionRate(), a_rlc->GetRecordSize(), a_rlc->GetSampleSize(),
                t_data_type_size, t_signed_vals,
                t_bit_depth, t_bits_right_aligned,
                &t_chan_vec );


        for( std::vector< unsigned >::const_iterator it = t_chan_vec.begin(); it != t_chan_vec.end(); ++it )
        {
            header->GetChannelHeaders()[ *it ].SetVoltageOffset( t_v_offset );
            header->GetChannelHeaders()[ *it ].SetVoltageRange( t_v_range );
//            header->GetChannelHeaders()[ *it ].SetDACGain( 1. );
            header->GetChannelHeaders()[ *it ].SetDACGain( t_dac_gain );
        }

        f_monarch->WriteHeader();

        f_stream = f_monarch->GetStream( t_stream_id );
        f_record = f_stream->GetStreamRecord();

        f_record_id = 0;
        f_record_time = 0;
        f_record_length = ( double )a_rlc->GetRecordSize() / ( 1.e-3 * a_rlc->GetAcquisitionRate() ); // in ns
        //f_record_n_bytes = n_channels * a_digitizer->DigitizerParams().data_type_size * a_rlc->GetRecordSize() * a_rlc->GetSampleSize();
        f_record_n_bytes = n_channels * t_data_type_size * a_rlc->GetRecordSize() * a_rlc->GetSampleSize();

        f_state = kPrepared;

        LINFO( lmclog, "Egg file <" << f_filename << "> is ready for records" );

        return true;
    }

    bool EggWriter::WriteRecord( const Signal* aSignal, bool is_new_acq )
    {

        if( f_state != kWriting && f_state != kPrepared )
        {
            LERROR( lmclog, "Egg file must be opened before writing records" );
            return false;
        }


        f_state = kWriting;

        LINFO( lmclog, "Writing record " << f_record_id );



        f_record->SetRecordId( f_record_id );

        f_record->SetTime( f_record_time );



        if( aSignal->GetState() != Signal::kDigital )
        {
            //LERROR( lmclog, "Signal is not digitized (state = " << aSignal->GetState() << "); no record was written" );
            //return false;
            ::memcpy( f_record->GetData(), (double*) aSignal->SignalTimeComplex(), f_record_n_bytes );
        } else {
 
			if( aSignal->GetDigitalIsSigned() ) ::memcpy( f_record->GetData(), reinterpret_cast< const monarch3::byte_type* >( aSignal->SignalDigitalS() ), f_record_n_bytes );
			else ::memcpy( f_record->GetData(), reinterpret_cast< const monarch3::byte_type* >( aSignal->SignalDigitalUS() ), f_record_n_bytes );

		}

        ++f_record_id;
        f_record_time += f_record_length;

        bool t_return = f_stream->WriteRecord( is_new_acq );
        //f_is_new_acq = false;
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

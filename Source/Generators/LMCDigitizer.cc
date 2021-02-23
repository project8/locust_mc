/*
 * LMCDigitizer.cc
 *
 *  Created on: Mar 3, 2014
 *      Author: nsoblath
 */

#include "LMCDigitizer.hh"

#include "logger.hh"


using scarab::get_calib_params;
using scarab::dig_calib_params;
using scarab::a2d;

namespace locust
{
    LOGGER( lmclog, "Digitizer" );

    MT_REGISTER_GENERATOR(Digitizer, "digitizer");

    Digitizer::Digitizer( const std::string& aName ) :
    		fRange( 2.e-8 ),
			fOffset( -1.e-8 ),
            Generator( aName ),
            fADCValuesSigned( false )
    {
        fRequiredSignalState = Signal::kTime;
		get_calib_params( 8, 1, -1.e-8, 2.e-8, false, &fParams );
    }

    Digitizer::~Digitizer()
    {
    }

    bool Digitizer::Configure( const scarab::param_node& aNode )
    {
        if( aNode.has( "adc-values-signed" ) )
            SetADCValuesSigned( aNode.get_value< bool >( "adc-values-signed", fADCValuesSigned ) );

        unsigned bitDepth = aNode.get_value( "bit-depth", fParams.bit_depth );
        unsigned dataTypeSize = aNode.get_value( "data-type-size", fParams.data_type_size );
        fRange = aNode.get_value( "v-range", fParams.v_range );
        fOffset = aNode.get_value( "v-offset", fParams.v_offset );
        

        get_calib_params( bitDepth, dataTypeSize, fOffset, fRange, false, &fParams );

        LDEBUG( lmclog, "Digitizer calibration parameters set" );

        return true;
    }


    void Digitizer::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

    const dig_calib_params& Digitizer::DigitizerParams() const
    {
        return fParams;
    }

    dig_calib_params& Digitizer::DigitizerParams()
    {
        return fParams;
    }

    bool Digitizer::DoGenerate( Signal* aSignal )
    {
  
        unsigned nchannels = fNChannels;
        unsigned signalSize = aSignal->TimeSize();
        unsigned signalSizeComplex = 2*aSignal->TimeSize()*nchannels;

        double* analogData = aSignal->SignalTime();
//        uint64_t* digitizedData = new uint64_t[ signalSize ];

        if( fADCValuesSigned )
        {

            int8_t* digitizedData = new int8_t[ signalSizeComplex ];

            for (unsigned ch = 0; ch < nchannels; ++ch)
            {
            for( unsigned index = 0; index < signalSize; ++index )
            {
                
                digitizedData[2*ch*signalSize + index*2 ] = a2d< double, int8_t >( aSignal->SignalTimeComplex()[ch*signalSize + index ][0], &fParams );
                digitizedData[2*ch*signalSize + index*2+1 ] = a2d< double, int8_t >( aSignal->SignalTimeComplex()[ch*signalSize + index ][1], &fParams );

                if ((int(digitizedData[ 2*ch*signalSize + index*2 ]) == 0 ) || ( int(digitizedData[ 2*ch*signalSize + index*2 ] == 255)))
                {
                    LERROR(lmclog,"Digitizer range limit.\n");
                    printf("Analog data at index %d channel %d is %g\n", index, ch, aSignal->SignalTimeComplex()[ch*signalSize + index ][0]);
                    printf("Digitized data at index %d channel %d is %d\n", index, ch, digitizedData[2*ch*signalSize + index*2 ]);
                	throw 1;
                	return false;
                }

                if( index < 10 )
                {

                    LWARN( lmclog, "digitizing channel " << ch << ": " << index << " I: " << aSignal->SignalTimeComplex()[ch*signalSize + index ][0] << " --> " << (int) digitizedData[2*ch*signalSize + index*2 ] );
                    LWARN( lmclog, "digitizing channel " << ch << ": " << index << " Q: " << aSignal->SignalTimeComplex()[ch*signalSize + index ][1] << " --> " << (int) digitizedData[2*ch*signalSize + index*2+1 ] );
                }

            } // signalsize
            } // channels
            aSignal->ToDigital( digitizedData, signalSizeComplex );
        }  // fADCValuesSigned
        else  // unsigned
        {
               
                uint8_t* digitizedData = new uint8_t[ signalSizeComplex ];

                for (unsigned ch = 0; ch < nchannels; ++ch)
                {
                for( unsigned index = 0; index < signalSize; ++index )
                {

                    digitizedData[2*ch*signalSize + index*2 ] = a2d< double, uint8_t >( aSignal->SignalTimeComplex()[ch*signalSize + index ][0], &fParams );
                    digitizedData[2*ch*signalSize + index*2+1 ] = a2d< double, uint8_t >( aSignal->SignalTimeComplex()[ch*signalSize + index ][1], &fParams );

                    if ((int(digitizedData[ 2*ch*signalSize + index*2 ]) == 0 ) || ( int(digitizedData[ 2*ch*signalSize + index*2 ] == 255)))
                    {
                        LERROR(lmclog,"Digitizer range limit.\n");
                        printf("Analog data at index %d channel %d is %g\n", index, ch, aSignal->SignalTimeComplex()[ch*signalSize + index ][0]);
                        printf("Digitized data at index %d channel %d is %d\n", index, ch, digitizedData[2*ch*signalSize + index*2 ]);
                    	throw 1;
                    	return false;
                    }

                    if( index < 20 )
                    {
                        LWARN( lmclog, "digitizing channel " << ch << ": " << index << " I: " << aSignal->SignalTimeComplex()[ch*signalSize + index ][0] << " --> " << (int) digitizedData[2*ch*signalSize + index*2 ] );
                        LWARN( lmclog, "digitizing channel " << ch << ": " << index << " Q: " << aSignal->SignalTimeComplex()[ch*signalSize + index ][1] << " --> " << (int) digitizedData[2*ch*signalSize + index*2+1 ] );
                    }

                } // signalsize
                } // channels
                aSignal->ToDigital( digitizedData, signalSizeComplex );

        }  // unsigned

        return true;

    }

} /* namespace locust */

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
            Generator( aName ),
            fADCValuesSigned( false )
    {
        fRequiredSignalState = Signal::kTime;

        get_calib_params( 8, 1, -3.e-6, 6.e-6, false, &fParams );  // if Gaussian noise is included.

//        get_calib_params( 8, 1, -1.e-7, 2.e-7, &fParams );  // if Gaussian noise is not included.
    }

    Digitizer::~Digitizer()
    {
    }

    bool Digitizer::Configure( const ParamNode* aNode )
    {
        if( aNode == NULL ) return true;

        unsigned bitDepth = aNode->GetValue( "bit-depth", fParams.bit_depth );
        unsigned dataTypeSize = aNode->GetValue( "data-type-size", fParams.data_type_size );
        double vRange = aNode->GetValue( "v-range", fParams.v_range );
        double vMin = aNode->GetValue( "v-offset", fParams.v_offset );

        get_calib_params( bitDepth, dataTypeSize, vMin, vRange, false, &fParams );

        DEBUG( lmclog, "Digitizer calibration parameters set" );

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

    bool Digitizer::DoGenerate( Signal* aSignal ) const
    {
        unsigned signalSize = aSignal->TimeSize();

        double* analogData = aSignal->SignalTime();
//        uint64_t* digitizedData = new uint64_t[ signalSize ];
        if( fADCValuesSigned )
        {
            int8_t* digitizedData = new int8_t[ signalSize ];

            for( unsigned index = 0; index < signalSize; ++index )
            {
                digitizedData[ index ] = a2d< double, int8_t >( analogData[ index ], &fParams );
                if( index < 100 )
                {
                    WARN( lmclog, "digitizing: " << index << ": " << analogData[ index ] << " --> " << (int) digitizedData[ index ] );  // pls added (int)
        //            printf("digitized data is %x\n", digitizedData[index]);
        //            getchar();
                }
            }
            aSignal->ToDigital( digitizedData, signalSize );
        }
        else
        {
            uint8_t* digitizedData = new uint8_t[ signalSize ];

            for( unsigned index = 0; index < signalSize; ++index )
            {
                digitizedData[ index ] = a2d< double, uint8_t >( analogData[ index ], &fParams );
                if( index < 100 )
                {
                    WARN( lmclog, "digitizing: " << index << ": " << analogData[ index ] << " --> " << (unsigned) digitizedData[ index ] );  // pls added (int)
        //            printf("digitized data is %x\n", digitizedData[index]);
        //            getchar();
                }
            }
            aSignal->ToDigital( digitizedData, signalSize );
        }

        return true;
    }

} /* namespace locust */

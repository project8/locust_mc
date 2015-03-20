/*
 * LMCDigitizer.cc
 *
 *  Created on: Mar 3, 2014
 *      Author: nsoblath
 */

#include "LMCDigitizer.hh"

#include "../Core/LMCLogger.hh"

namespace locust
{
    LMCLOGGER( lmclog, "Digitizer" );

    MT_REGISTER_GENERATOR(Digitizer, "digitizer");

    Digitizer::Digitizer( const std::string& aName ) :
            Generator( aName )
    {
        fRequiredSignalState = Signal::kTime;

        get_calib_params( 8, 1, -1.e-7, 2.e-7, &fParams );
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
        double vMin = aNode->GetValue( "v-min", fParams.v_min );

        get_calib_params( bitDepth, dataTypeSize, vMin, vRange, &fParams );

        LMCDEBUG( lmclog, "Digitizer calibration parameters set:\n" << fParams );

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
        uint8_t* digitizedData = new uint8_t[ signalSize ];

        for( unsigned index = 0; index < signalSize; ++index )
        {
            digitizedData[ index ] = a2d( analogData[ index ], &fParams );
            if( index < 100 )
            {
                LMCWARN( lmclog, "digitizing: " << index << ": " << analogData[ index ] << " --> " << digitizedData[ index ] );
            }
        }

        aSignal->ToDigital( digitizedData, signalSize );
        return true;
    }

} /* namespace locust */

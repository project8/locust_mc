/*
 * LMCDigitizer.cc
 *
 *  Created on: Mar 3, 2014
 *      Author: nsoblath
 */

#include "LMCDigitizer.hh"

namespace locust
{
    MT_REGISTER_GENERATOR(Digitizer, "digitizer");

    Digitizer::Digitizer( const std::string& aName ) :
            Generator( aName )
    {
        get_calib_params( 8, 1, -0.25, 0.5, &fParams );
    }

    Digitizer::~Digitizer()
    {
    }

    void Digitizer::Configure( const ParamNode* aNode )
    {
        unsigned bitDepth = aNode->GetValue( "bit-depth", fParams.bit_depth );
        unsigned dataTypeSize = aNode->GetValue( "data-type-size", fParams.data_type_size );
        double vRange = aNode->GetValue( "v-range", fParams.v_range );
        double vMin = aNode->GetValue( "v-min", fParams.v_min );

        get_calib_params( bitDepth, dataTypeSize, vMin, vRange, &fParams );

        return;
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

    void Digitizer::Generate( Signal* aSignal ) const
    {
        unsigned signalSize = aSignal->TimeSize();

        double* analogData = aSignal->SignalTime();
        uint64_t* digitizedData = new uint64_t[ signalSize ];

        for( unsigned index = 0; index < signalSize; ++index )
        {
            digitizedData[ index ] = a2d( analogData[ index ], &fParams );
        }

        aSignal->ToDigital( digitizedData, signalSize );
        return;
    }

} /* namespace locust */

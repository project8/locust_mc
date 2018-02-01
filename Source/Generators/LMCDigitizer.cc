/*
 * LMCDigitizer.cc
 *
 *  Created on: Mar 3, 2014
 *      Author: nsoblath
 */

#include "LMCDigitizer.hh"

#include "logger.hh"
#include "LMCSimulationController.hh"


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

	//					        get_calib_params( 8, 1, -3.e-6, 6.e-6, false, &fParams );  // if Gaussian noise is included.

						      get_calib_params( 8, 1, -1.e-7, 2.e-7, false, &fParams );  // if Gaussian noise is not included.
    }

    Digitizer::~Digitizer()
    {
    }

    bool Digitizer::Configure( const scarab::param_node* aNode )
    {
        if( aNode == NULL ) return true;

        unsigned bitDepth = aNode->get_value( "bit-depth", fParams.bit_depth );
        unsigned dataTypeSize = aNode->get_value( "data-type-size", fParams.data_type_size );
        double vRange = aNode->get_value( "v-range", fParams.v_range );
        double vMin = aNode->get_value( "v-offset", fParams.v_offset );

        get_calib_params( bitDepth, dataTypeSize, vMin, vRange, false, &fParams );

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

    bool Digitizer::DoGenerate( Signal* aSignal ) const
    {

    	bool IQStream = true;  // this could eventually be a parameter in the json file.

    	SimulationController SimulationController1;
        const unsigned nchannels = SimulationController1.GetNChannels();


        unsigned signalSize = aSignal->TimeSize();
        unsigned signalSizeComplex = 2*aSignal->TimeSize()*nchannels;

        double* analogData = aSignal->SignalTime();
//        uint64_t* digitizedData = new uint64_t[ signalSize ];

/*
        std::complex<double>* analogDataComplex;
        analogDataComplex = (std::complex<double> *) malloc(nchannels*signalSize * sizeof(std::complex<double>));
        memcpy( analogDataComplex, aSignal->SignalTimeComplex(), nchannels*signalSize * sizeof( std::complex<double> ) );
*/

        if( fADCValuesSigned )
        {
        if (!IQStream)
        {
            int8_t* digitizedData = new int8_t[ signalSize ];

            for( unsigned index = 0; index < signalSize; ++index )
            {
                digitizedData[ index ] = a2d< double, int8_t >( analogData[ index ], &fParams );
                if( index < 10 )
                {
                    LWARN( lmclog, "digitizing: " << index << ": " << analogData[ index ] << " --> " << (int) digitizedData[ index ] );  // pls added (int)
                }
            }
            aSignal->ToDigital( digitizedData, signalSize );
        }  // !IQStream
        else
        {
            int8_t* digitizedData = new int8_t[ signalSizeComplex ];
//            printf("signalSizeComplex is %d\n", signalSizeComplex); getchar();

            for (unsigned ch = 0; ch < nchannels; ++ch)
            {
            for( unsigned index = 0; index < signalSize; ++index )
            {
//            	printf("2*ch*signalSize+index*2 is %d\n", 2*ch*signalSize+index*2);
//                digitizedData[2*ch*signalSize + index*2 ] = a2d< double, int8_t >( analogDataComplex[ch*signalSize + index ].real(), &fParams );
//                digitizedData[2*ch*signalSize + index*2+1 ] = a2d< double, int8_t >( analogDataComplex[ch*signalSize + index ].imag(), &fParams );
                digitizedData[2*ch*signalSize + index*2 ] = a2d< double, int8_t >( aSignal->SignalTimeComplex()[ch*signalSize + index ][0], &fParams );
                digitizedData[2*ch*signalSize + index*2+1 ] = a2d< double, int8_t >( aSignal->SignalTimeComplex()[ch*signalSize + index ][1], &fParams );

                if( index < 10 )
                {

                    LWARN( lmclog, "digitizing channel " << ch << ": " << index << " I: " << aSignal->SignalTimeComplex()[ch*signalSize + index ][0] << " --> " << (int) digitizedData[2*ch*signalSize + index*2 ] );  // pls added (int)
                    LWARN( lmclog, "digitizing channel " << ch << ": " << index << " Q: " << aSignal->SignalTimeComplex()[ch*signalSize + index ][1] << " --> " << (int) digitizedData[2*ch*signalSize + index*2+1 ] );  // pls added (int)
                }
            }
            }
            aSignal->ToDigital( digitizedData, signalSizeComplex );
        }  // IQStream
        }
        else
        {
        	if (!IQStream)
        	{
            uint8_t* digitizedData = new uint8_t[ signalSize ];

            for( unsigned index = 0; index < signalSize; ++index )
            {
                digitizedData[ index ] = a2d< double, uint8_t >( analogData[ index ], &fParams );
                if( index < 10 )
                {
                    LWARN( lmclog, "digitizing: " << index << ": " << analogData[ index ] << " --> " << (unsigned) digitizedData[ index ] );  // pls added (int)
        //            printf("digitized data is %x\n", digitizedData[index]);
        //            getchar();
                }
            }
            aSignal->ToDigital( digitizedData, signalSize );
        	} // !IQStream
        	else
        	{
                uint8_t* digitizedData = new uint8_t[ signalSizeComplex ];
//                printf("signalSizeComplex is %d\n", signalSizeComplex); getchar();

                for (unsigned ch = 0; ch < nchannels; ++ch)
                {
                for( unsigned index = 0; index < signalSize; ++index )
                {
//                	printf("ch is %d, index is %d, 2*ch*signalSize+index*2 is %d\n", ch, index, 2*ch*signalSize+index*2);
//                    digitizedData[2*ch*signalSize + index*2 ] = a2d< double, uint8_t >( analogDataComplex[ch*signalSize + index ].real(), &fParams );
//                    digitizedData[2*ch*signalSize + index*2+1 ] = a2d< double, uint8_t >( analogDataComplex[ch*signalSize + index ].imag(), &fParams );
                    digitizedData[2*ch*signalSize + index*2 ] = a2d< double, uint8_t >( aSignal->SignalTimeComplex()[ch*signalSize + index ][0], &fParams );
                    digitizedData[2*ch*signalSize + index*2+1 ] = a2d< double, uint8_t >( aSignal->SignalTimeComplex()[ch*signalSize + index ][1], &fParams );

                	// fake data for debugging.
//                    digitizedData[2*ch*signalSize + index*2 ] = a2d< double, uint8_t >( 5.e-8, &fParams );
//                    digitizedData[2*ch*signalSize + index*2+1 ] = a2d< double, uint8_t >( 5.e-8, &fParams );


                    if( index < 20 )
                    {
                        LWARN( lmclog, "digitizing channel " << ch << ": " << index << " I: " << aSignal->SignalTimeComplex()[ch*signalSize + index ][0] << " --> " << (int) digitizedData[2*ch*signalSize + index*2 ] );  // pls added (int)
                        LWARN( lmclog, "digitizing channel " << ch << ": " << index << " Q: " << aSignal->SignalTimeComplex()[ch*signalSize + index ][1] << " --> " << (int) digitizedData[2*ch*signalSize + index*2+1 ] );  // pls added (int)
                    }
                }
                }
                aSignal->ToDigital( digitizedData, signalSizeComplex );

        	}  // IQStream
        }

        return true;
    }

} /* namespace locust */

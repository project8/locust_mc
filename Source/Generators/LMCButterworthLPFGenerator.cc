/*
 * LMCButterworthLPFGenerator.cc
 *
 *  Created on: Sept 9, 2016
 *      Author: plslocum after nsoblath
 */

#include "LMCButterworthLPFGenerator.hh"

#include "LMCGlobalsDeclaration.hh"
#include <deque>

#include "logger.hh"

using std::string;

namespace locust
{
    LOGGER( lmclog, "ButterworthLPFGenerator" );

    MT_REGISTER_GENERATOR(ButterworthLPFGenerator, "butterworth-lpf");

    ButterworthLPFGenerator::ButterworthLPFGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &ButterworthLPFGenerator::DoGenerateTime )
    {
        fRequiredSignalState = Signal::kTime;
    }

    ButterworthLPFGenerator::~ButterworthLPFGenerator()
    {
    }

    bool ButterworthLPFGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;
        return true;
    }

    void ButterworthLPFGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }


    bool ButterworthLPFGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }

    bool ButterworthLPFGenerator::DoGenerateTime( Signal* aSignal ) const
    {
    	// 8th order Butterworth filter with wc = 70.e6 Hz fs=200 MHz.

    	double* ai = new double[20];
    	double* bi = new double[20];
    	int N = 8;  // 8th order digital filter.
    	int M = 8;

        //Denominator Coefficients
        //ai[0]=1.0000;
        //ai[1]=-7.436415727064;
        //ai[2]=24.212528496827;
        //ai[3]=-45.08206974092;
        //ai[4]=52.5005503707;
        //ai[5]=-39.15720043493;
        //ai[6]=18.265819680576;
        //ai[7]=-4.872150915425;
        //ai[8]=0.568938280997;
        //ai[9]=0.0000000;

        ai[0]=1.0000;
        ai[1]=-6.8730110236;
        ai[2]=20.7365964965756;
        ai[3]=-35.864676348246;
        ai[4]=38.884394537568;
        ai[5]=-27.057833592167;
        ai[6]=11.799329601546;
        ai[7]=-2.947769270203;
        ai[8]=0.322972809305;
        ai[9]=0.0000000;

        //Numerator coefficients
        //bi[0]=6.3465677158e-11;
        //bi[1]=5.0772541726e-10;
        //bi[2]=1.77703896042e-9;
        //bi[3]=3.554077920853e-9;
        //bi[4]=4.442597401066e-9;
        //bi[5]=3.554077920853e-9;
        //bi[6]=1.77703896042e-9;
        //bi[7]=5.0772541726e-10;
        //bi[8]=6.3465677158e-11;
        //bi[9]=0.0000000;
        
        bi[0]=1.25420323016e-8;
        bi[1]=1.00336258412e-7;
        bi[2]=3.511769044450e-7;
        bi[3]=7.0235380889e-7;
        bi[4]=8.77942261112e-7;
        bi[5]=7.0235380889e-7;
        bi[6]=3.511769044450e-7;
        bi[7]=1.00336258412e-7;
        bi[8]=1.25420323016e-8;
        bi[9]=0.0000000;


        //Toss out useless normalization factor
        for(unsigned int index = 0 ; index < M ; ++index)
        {
            ai[index]=ai[index+1];
        }

        double *FilteredSignal = new double[aSignal->TimeSize()]; // output signal.
        std::deque<double> InternalSignal(N); // temporary signal

        //Initialize
        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
            FilteredSignal[index] = 0.;
        
        for( unsigned index = 0; index < N; ++index )
            InternalSignal[index]=0.0;

        for( unsigned index = 0; index < aSignal->TimeSize() ; ++index)
        {
            FilteredSignal[index]=bi[0]*aLongSignal[index]+InternalSignal[0];
            InternalSignal.pop_front();
            InternalSignal.push_back(0.);
            
            for(unsigned j = 0; j < N; j++)
            {
                InternalSignal[j]=InternalSignal[j]+bi[j+1]*aLongSignal[index];
            }

            for(unsigned j = 0; j < M; j++)
            {
                InternalSignal[j]=InternalSignal[j]-ai[j]*FilteredSignal[index];
            }


        }

        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
            aLongSignal[index] = FilteredSignal[index];
        }

        delete FilteredSignal;
        return true;
    }

    bool ButterworthLPFGenerator::DoGenerateFreq( Signal* aSignal ) const
    {

        return true;
    }

} /* namespace locust */

/*
 * LMCVisitor.hh
 *
 *  Created on: Feb 7, 2014
 *      Author: nsoblath
 */

#ifndef LMCVISITOR_HH_
#define LMCVISITOR_HH_

namespace locust
{
    class Generator;
    class GaussianNoiseGenerator;
    class TestSignalGenerator;
    class BasebandTrackGenerator;
    class LowPassFilterFFTGenerator;
    class Digitizer;
    class TrappedElectronGenerator;
    class KassSignalGenerator;
    class DecimateSignalGenerator;

    class GeneratorVisitor
    {
        public:
            GeneratorVisitor();
            virtual ~GeneratorVisitor();

            virtual void Visit( const Generator* );

            virtual void Visit( const KassSignalGenerator* ) = 0;
            virtual void Visit( const GaussianNoiseGenerator* ) = 0;
            virtual void Visit( const TrappedElectronGenerator* ) = 0;
            virtual void Visit( const TestSignalGenerator* ) = 0;
            virtual void Visit( const BasebandTrackGenerator* ) = 0;
            virtual void Visit( const LowPassFilterFFTGenerator* ) = 0;
            virtual void Visit( const DecimateSignalGenerator* ) = 0;
            virtual void Visit( const Digitizer* ) = 0;
    };


} /* namespace locust */

#endif /* LMCVISITOR_HH_ */

/*
 * LMCRunLengthCalculator.cc
 *
 *  Created on: Feb 7, 2014
 *      Author: nsoblath
 */

#include "LMCRunLengthCalculator.hh"

#include "LMCGenerator.hh"
#include "logger.hh"
#include "LMCParam.hh"

namespace locust
{
    LOGGER( lmclog, "RunLengthCalculator" );

    RunLengthCalculator::RunLengthCalculator() :
            GeneratorVisitor(),
            fState( kUnknown ),
            fNRecords( 1 ),
            fDuration( 1. ),
            fFirstGenerator( NULL ),
            fByGeneratorsState( kUnknown ),
            fByGeneratorsNRecords( 1 ),
            fByGeneratorsDuration( 1. ),
            fRecordSize( 4194304 ),
            fBinWidth( 5.e-9 ),
            fAcquisitionRate( 200. )
    {
    }

    RunLengthCalculator::~RunLengthCalculator()
    {
    }

    bool RunLengthCalculator::Configure( const ParamNode* aNode )
    {
        if( aNode == NULL ) return false;

        // first, configure items that will affect the method for calculating the run length

        if( aNode->Has( "n-records" ) )
            SetNRecords( aNode->GetValue< unsigned >( "n-records" ) );

        if( aNode->Has( "duration" ) )
            SetDuration( aNode->GetValue< double >( "duration" ) );

        // next, configure items that are needed no matter what.

        SetRecordSize( aNode->GetValue< unsigned >( "record-size", fRecordSize ) );

        if( aNode->Has( "acquisition-rate" ) )
            SetAcquisitionRate( aNode->GetValue< double >( "acquisition-rate" ) );
        if( aNode->Has( "bin-width" ) )
            SetBinWidth( aNode->GetValue< double >( "bin-width" ) );

        return true;
    }

    bool RunLengthCalculator::VisitGenerators()
    {
        if( fFirstGenerator == NULL )
        {
            WARN( lmclog, "First generator has not been set" );
            return false;
        }

        const Generator* nextGenerator = fFirstGenerator;
        while( nextGenerator != NULL )
        {
            nextGenerator->Accept( this );
            nextGenerator = nextGenerator->GetNextGenerator();
        }

        return true;
    }

    void RunLengthCalculator::Visit( const KassSignalGenerator* )
    {
        // nothing to see here, move along, please
        return;
    }
    void RunLengthCalculator::Visit( const FreeFieldSignalGenerator* )
    {
        // nothing to see here, move along, please
        return;
    }

    void RunLengthCalculator::Visit( const TrappedElectronGenerator* )
    {
        // nothing to see here, move along, please
        return;
    }

    void RunLengthCalculator::Visit( const LowPassFilterFFTGenerator* )
    {
        // nothing to see here, move along, please
        return;
    }

    void RunLengthCalculator::Visit( const GaussianNoiseGenerator* )
    {
        // nothing to see here, move along, please
        return;
    }

    void RunLengthCalculator::Visit( const Digitizer* )
    {
        // nothing to see here, move along, please
        return;
    }
    void RunLengthCalculator::Visit( const DecimateSignalGenerator* )
     {
         // nothing to see here, move along, please
         return;
     }



    void RunLengthCalculator::Visit( const TestSignalGenerator* )
    {
        // nothing to see here, move along, please
        return;
    }

  void RunLengthCalculator::Visit( const BasebandTrackGenerator* )
    {
        // nothing to see here, move along, please
        return;
    }

    bool RunLengthCalculator::CalculateRunLength()
    {
        switch( fState )
        {
            case kDecided:
                WARN( lmclog, "Run length already calculated; use Reset function to recalculate" );
                return false;
                break;
            case kByRecords:
                fDuration = fNRecords * fBinWidth * fRecordSize;
                fState = kDecided;
                INFO( lmclog, "Run length determined by the number of records" );
                break;
            case kByDuration:
                fNRecords = (unsigned)( fDuration / double(fBinWidth * fRecordSize) );
                fDuration = fNRecords * fBinWidth * fRecordSize;
                fState = kDecided;
                INFO( lmclog, "Run length determined by the duration, rounded down to an integer number of records" );
                break;
            case kByGenerators:
                if( fByGeneratorsState == kByRecords )
                {
                    fNRecords = fByGeneratorsNRecords;
                    fDuration = fNRecords * fBinWidth * fRecordSize;
                    fState = kDecided;
                    INFO( lmclog, "Run length determined by generators setting the number of records" );
                    break;
                }
                else if( fByGeneratorsState == kByDuration )
                {
                    fDuration = fByGeneratorsDuration;
                    fNRecords = (unsigned)( fDuration / double(fBinWidth * fRecordSize) );
                    fDuration = fNRecords * fBinWidth * fRecordSize;
                    fState = kDecided;
                    INFO( lmclog, "Run length determined by generators setting the duration" );
                    break;
                }
                WARN( lmclog, "Unable to set parameters by generators (by-generators state: " << fByGeneratorsState << ")" );
                return false;
                break;
            case kUnknown:
                WARN( lmclog, "Run length parameters have not been supplied" );
                return false;
                break;
            default:
                WARN( lmclog, "Unrecognized state: " << fState );
                return false;
                break;
        }

        INFO( lmclog, "Final run length parameters:\n" <<
                 "\t\tNumber of records: " << fNRecords << '\n' <<
                 "\t\tDuration: " << fDuration << " seconds\n" <<
                 "\t\tRecord size: " << fRecordSize << " samples\n" <<
                 "\t\tAcquisition rate: " << fAcquisitionRate << " MHz\n" <<
                 "\t\tBin width: " << fBinWidth << " s");
        return true;
    }

    RunLengthCalculator::RunLengthState RunLengthCalculator::GetState() const
    {
        return fState;
    }

    void RunLengthCalculator::SetRunLengthState( RunLengthState state )
    {
        if( state < fState ) fState = state;
        return;
    }

    void RunLengthCalculator::OverrideRunLengthState( RunLengthState state )
    {
        fState = state;
        return;
    }

    unsigned RunLengthCalculator::GetNRecords() const
    {
        return fNRecords;
    }

    void RunLengthCalculator::SetNRecords( unsigned recs )
    {
        fNRecords = recs;
        SetRunLengthState( kByRecords );
        return;
    }

    double RunLengthCalculator::GetDuration() const
    {
        return fDuration;
    }

    void RunLengthCalculator::SetDuration( double duration )
    {
        fDuration = duration;
        SetRunLengthState( kByDuration );
        return;
    }

    void RunLengthCalculator::SetFirstGenerator( const Generator* firstGen )
    {
        fFirstGenerator = firstGen;
        return;
    }

    RunLengthCalculator::RunLengthState RunLengthCalculator::GetByGeneratorsState() const
    {
        return fByGeneratorsState;
    }

    unsigned RunLengthCalculator::GetByGeneratorsNRecords() const
    {
        return fByGeneratorsNRecords;
    }

    double RunLengthCalculator::GetByGeneratorsDuration() const
    {
        return fByGeneratorsDuration;
    }

    unsigned RunLengthCalculator::GetRecordSize() const
    {
        return fRecordSize;
    }

    void RunLengthCalculator::SetRecordSize( unsigned size )
    {
        fRecordSize = size;
        return;
    }

    double RunLengthCalculator::GetAcquisitionRate() const
    {
        return fAcquisitionRate;
    }

    void RunLengthCalculator::SetAcquisitionRate( double ar )
    {
        fAcquisitionRate = ar;
        fBinWidth = 1. / (ar * 1.e6);
        return;
    }

    double RunLengthCalculator::GetBinWidth() const
    {
        return fBinWidth;
    }

    void RunLengthCalculator::SetBinWidth( double bw )
    {
        fBinWidth = bw;
        fAcquisitionRate = 1. / (bw * 1.e6);
        return;
    }



} /* namespace locust */

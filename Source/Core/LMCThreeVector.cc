#include "LMCThreeVector.hh"

namespace locust
{

    const LMCThreeVector LMCThreeVector::sInvalid( NAN, NAN, NAN );
    const LMCThreeVector LMCThreeVector::sZero( 0., 0., 0. );

    const LMCThreeVector LMCThreeVector::sXUnit( 1., 0., 0. );
    const LMCThreeVector LMCThreeVector::sYUnit( 0., 1., 0. );
    const LMCThreeVector LMCThreeVector::sZUnit( 0., 0., 1. );

    LMCThreeVector::LMCThreeVector()
    {
        fData[0] = 0.;
        fData[1] = 0.;
        fData[2] = 0.;
    }
    //LMCThreeVector::~LMCThreeVector()
    //{
    //}

}

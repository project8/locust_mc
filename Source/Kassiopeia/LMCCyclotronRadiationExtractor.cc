/*
 * LMCCyclotronRadiationExtractor.cc
 *
 *  Created on: Mar 13, 2016
 *      Author: nsoblath
 */

#include "LMCCyclotronRadiationExtractor.hh"

using namespace Kassiopeia;

namespace locust
{

    CyclotronRadiationExtractor::CyclotronRadiationExtractor()
    {
    }

    CyclotronRadiationExtractor::CyclotronRadiationExtractor( const CyclotronRadiationExtractor& aOrig )
    {
    }

    CyclotronRadiationExtractor::~CyclotronRadiationExtractor()
    {
    }

    CyclotronRadiationExtractor* CyclotronRadiationExtractor::Clone() const
    {
        return new CyclotronRadiationExtractor( *this );
    }


    bool CyclotronRadiationExtractor::ExecutePreStepModification( KSParticle& anInitialParticle, KSParticleQueue& aQueue )
    {
        return true;
    }

    bool CyclotronRadiationExtractor::ExecutePostStepModifcation( KSParticle& anInitialParticle, KSParticle& aFinalParticle, KSParticleQueue& aQueue )
    {
        return true;
    }


    void CyclotronRadiationExtractor::InitializeComponent()
    {
    }

    void CyclotronRadiationExtractor::DeinitializeComponent()
    {
    }

    void CyclotronRadiationExtractor::PullDeupdateComponent()
    {
    }
    void CyclotronRadiationExtractor::PushDeupdateComponent()
    {
    }


} /* namespace locust */

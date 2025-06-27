
#include "LMCPowerNormFieldCalculator.hh"
#include "logger.hh"
#include "KSParticleFactory.h"
#include <algorithm>

#include "KSInteractionsMessage.h"
#include <limits>

using std::numeric_limits;

using namespace Kassiopeia;
namespace locust
{

    LOGGER( lmclog, "PowerNormFieldCalculator" );

    PowerNormFieldCalculator::PowerNormFieldCalculator() :
        fNFilterBinsRequired( 0 ),
        fTFReceiverHandler( NULL ),
        fAnalyticResponseFunction( 0 ),
        fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }
    PowerNormFieldCalculator::PowerNormFieldCalculator( const PowerNormFieldCalculator& aCopy ) :
        fNFilterBinsRequired( 0 ),
        fTFReceiverHandler( NULL ),
        fAnalyticResponseFunction( 0 ),
        fInterface( aCopy.fInterface )
    {
    }
    PowerNormFieldCalculator* PowerNormFieldCalculator::Clone() const
    {
        return new PowerNormFieldCalculator( *this );
    }
    PowerNormFieldCalculator::~PowerNormFieldCalculator()
    {
        if (fTFReceiverHandler != NULL)
        {
            delete fTFReceiverHandler;
        }
        if (fAnalyticResponseFunction != NULL)
        {
            delete fAnalyticResponseFunction;
        }
    }


}

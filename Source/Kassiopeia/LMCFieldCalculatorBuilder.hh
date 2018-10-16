/*
 * LMCFieldCalculatorBuilder.hh
 *
 *  Created on: Oct 4, 2018
 *      Author: pslocum
 */

#ifndef LOCUST_LMCFIELDCALCULATORBUILDER_HH_
#define LOCUST_LMCFIELDCALCULATORBUILDER_HH_

#include "KComplexElement.hh"
#include "LMCFieldCalculator.hh"


using namespace locust;
namespace katrin
{
    /*!
     @class 
     @author P. L. Slocum

     @brief bindings for FieldCalculator for the locust xml file

     @details

     Available configuration options:
     - "name": string -- name of xml object
     - "P8Phase": int -- phase of P8 experiment. Currently 1, 2, 3 are supported
    */

    typedef KComplexElement< locust::FieldCalculator > FieldCalculatorBuilder;

    template< >
    inline bool FieldCalculatorBuilder::AddAttribute(KContainer *aContainer)
    {
        if( aContainer->GetName() == "name" )
        {
            aContainer->CopyTo( fObject, &KNamed::SetName );
            return true;
        }

       
        return false;
    }

} /* namespace katrin */

#endif /* LOCUST_LMCFIELDCALCULATORBUILDER_HH_ */


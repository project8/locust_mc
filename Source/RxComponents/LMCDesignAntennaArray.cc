/*
 * LMCDesignAntennaArray.cc
 *
 *  Created on: Oct 29, 2020
 *      Author: Arina Telles
 */

#include "LMCDesignAntennaArray.hh"
using std::string;


namespace locust
{

    LOGGER( lmclog, "DesignAntennaArray" );

    DesignAntennaArray::DesignAntennaArray():
            fImpedanceTransformation( 0.358 ) // sqrt(50./390.), where 390 ohms is the 5-slot TF input array impedance
    {
    }

    DesignAntennaArray::~DesignAntennaArray()
    {
    }

    bool DesignAntennaArray::Configure( const scarab::param_node& aParam )
    {

        if( !PowerCombiner::Configure(aParam))
        {
            LERROR(lmclog,"Error configuring PowerCombiner class from DesignAntennaArray subclass");
            return false;
        }

        if( aParam.has( "impedance-transformation" ) )
        {
            fImpedanceTransformation = aParam["impedance-transformation"]().as_double();
        }

        return true;
    }


    Receiver* DesignAntennaArray::ChooseElement()
    {
        DesignAntenna* aDesignAntenna = new DesignAntenna;
        return aDesignAntenna;
    }



} /* namespace locust */

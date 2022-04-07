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
            fImpedanceTransformation( 1.0 ) // default to 50 ohms
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

        if( aParam.has( "theta-pattern-exponent" ) )
        {
            fThetaPatternParameter = aParam["theta-pattern-exponent"]().as_double();
        }

        if( aParam.has( "phi-pattern-exponent" ) )
        {
            fPhiPatternParameter = aParam["phi-pattern-exponent"]().as_double();
        }


        SetVoltageDampingFactors();

        return true;
    }

    bool DesignAntennaArray::SetVoltageDampingFactors()
    {

        for (unsigned z_index=0; z_index<GetNElementsPerStrip(); z_index++)
        {
            double aFactor = fImpedanceTransformation;
            // no scaling by the number of elements here
            // transfer function should be made for 1 element.
            SetDampingFactor(z_index, aFactor);
        }

        return true;
    } 


    Receiver* DesignAntennaArray::ChooseElement()
    {
        DesignAntenna* aDesignAntenna = new DesignAntenna;
        aDesignAntenna->SetFieldPatternExponents(fThetaPatternParameter, fPhiPatternParameter);
        return aDesignAntenna;
    }



} /* namespace locust */

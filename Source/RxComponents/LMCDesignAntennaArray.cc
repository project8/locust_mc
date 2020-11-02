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
            //Adhoc scaling in case the slotted waveguide has fewer/more slots than 5
            aFactor = aFactor/sqrt(GetNElementsPerStrip()/5.);
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

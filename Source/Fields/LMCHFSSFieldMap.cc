/*
 * LMCHFSSFieldMap.cc
 *
 *  Created on: Jun 9, 2021
 *      Author: pslocum
 */

#include "LMCHFSSFieldMap.hh"

namespace locust
{
    LOGGER( lmclog, "HFSSFieldMap" );
    HFSSFieldMap::HFSSFieldMap():
    fInterface( KLInterfaceBootstrapper::get_instance()->GetInterface() )
    {
    }
    HFSSFieldMap::~HFSSFieldMap() {}

    void HFSSFieldMap::ReadHFSSFile()
    {
    	return;
    }


    double HFSSFieldMap::Integrate(int l, int m, int n, bool teMode, bool eField)
    {

    	std::vector<double> aField;
    	double xPozar, yPozar, xKass, yKass = 0.;
    	double dX = fInterface->fX/GetNPixels();
    	double dY = fInterface->fY/GetNPixels();
    	double tArea = 0.;
    	double tIntegral = 0.;

    	for (unsigned i=0; i<GetNPixels(); i++)
    		for (unsigned j=0; j<GetNPixels(); j++)
    		{
    	    	xPozar = (double)i*dX;
    	    	yPozar = (double)j*dY;
   	    		xKass = xPozar - fInterface->fX/2.;
   	    		yKass = yPozar - fInterface->fY/2.;

    	    	if (teMode)
    	    	{
    	    		if (eField)
    	    		{
    	    		    aField = TE_E(m, n, xKass, yKass);
    	    		}
    	    		else
    	    		{
    	    			aField = TE_H(m, n, xKass, yKass);
    	    		}
    	    	}
    	    	else
    	    	{
    	    		if (eField)
    	    		{
    	    			aField = TM_E(m, n, xKass, yKass);
    	    		}
    	    		else
    	    		{
    	    			aField = TM_H(m, n, xKass, yKass);
    	    		}
    	    	}

    	    	double aFieldMagSq = 0.;
    	    	auto it = aField.begin();
    	    	while (it != aField.end())
    	    	{
		    		if (!isnan(*it))
		    			aFieldMagSq += (*it)*(*it);
    	    		*it++;
    	    	}

    			tIntegral += aFieldMagSq*dX*dY;
//    		    tArea += dX*dY;  // sanity check area integral.
    		}
//    	printf("tArea is %g\n", tArea); getchar();
    	return tIntegral;
    }

    std::vector<double> HFSSFieldMap::TE_E(int m, int n, double xKass, double yKass) const
    {

    	double x = xKass + fInterface->fX/2.;
    	double y = yKass + fInterface->fY/2.;

    	// from Pozar
    	std::vector<double> TE_E;
        return TE_E;
    }

    std::vector<double> HFSSFieldMap::TE_H(int m, int n, double xKass, double yKass) const
    {
    	double x = xKass + fInterface->fX/2.;
    	double y = yKass + fInterface->fY/2.;

    	// from Pozar
    	std::vector<double> TE_H;
        return TE_H;
    }

    std::vector<double> HFSSFieldMap::TM_E(int m, int n, double xKass, double yKass) const
    {
    	double x = xKass + fInterface->fX/2.;
    	double y = yKass + fInterface->fY/2.;

    	// from Pozar
    	std::vector<double> TM_E;
        return TM_E;
    }

    std::vector<double> HFSSFieldMap::TM_H(int m, int n, double xKass, double yKass) const
    {
    	double x = xKass + fInterface->fX/2.;
    	double y = yKass + fInterface->fY/2.;

    	// from Pozar
    	std::vector<double> TM_H;
        return TM_H;
    }




} /* namespace locust */


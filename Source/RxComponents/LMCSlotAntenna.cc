/*
 * LMCSlotAntenna.cc
 *
 *  Created on: Dec 19, 2019
 *      Author: pslocum
 */

#include "LMCSlotAntenna.hh"


namespace locust
{

    SlotAntenna::SlotAntenna()
    {
    }

    SlotAntenna::~SlotAntenna()
    {
    }

    void SlotAntenna::RxSayHello()
    {
      	printf("slot says hello\n");
      	getchar();
    }

    double SlotAntenna::GetPatternFactor(LMCThreeVector incidentKVector, Receiver currentElement)
    {
    	// This is the dipole aoi factor only.  It is not the dot product with the co-pol direction.

    	double incidentNormal = (-1.0) * incidentKVector.Dot(currentElement.GetNormalDirection());
    	double incidentCoPol = currentElement.GetPolarizationDirection().Dot(incidentKVector);
    	double incidentCrossPol = currentElement.GetCrossPolarizationDirection().Dot(incidentKVector);
    	double tTheta = LMCConst::Pi()/2.;
        if (incidentCoPol>0.) tTheta = atan(incidentNormal/fabs(incidentCoPol));
        if (incidentCoPol<0.) tTheta = LMCConst::Pi() - atan(incidentNormal/fabs(incidentCoPol));

        // dipole theta dependence, normalized to 1.0 for normal incidence
        double dipoleThetaFactor = 0.0;
        if (fabs(tTheta)>0.) dipoleThetaFactor = cos((LMCConst::Pi()/2.)*cos(tTheta))/sin(tTheta);

        double tPhi = 0.;
    	if (fabs(incidentNormal)>0.) tPhi = atan(incidentCrossPol/incidentNormal);

    	// dipole donut pinch from HFSS, normalized to 1.0 for normal incidence.
    	double dipolePhiPinchFactor = cos(tPhi);

    	return dipoleThetaFactor * dipolePhiPinchFactor;
    }


} /* namespace locust */

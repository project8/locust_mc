/*
 * LMCDesignAntenna.cc
 *
 *  Created on: Oct 26, 2020
 *      Author: Arina Telles
 */

#include "LMCDesignAntenna.hh"


namespace locust
{

    DesignAntenna::DesignAntenna()
    {
    }

    DesignAntenna::~DesignAntenna()
    {
    }

    void DesignAntenna::RxSayHello()
    {
      	printf("Design Antenna says hello\n");
      	getchar();
    }

    double DesignAntenna::GetPatternFactor(LMCThreeVector incidentKVector, Receiver currentElement)
    {
    	// This is the aoi factor only.  It is not the dot product with the co-pol direction.

        // Calculate components of incident radiation along normal, co-pol, and cross-pol directions
    	double incidentNormal = (-1.0) * incidentKVector.Dot(currentElement.GetNormalDirection());
    	double incidentCoPol = currentElement.GetPolarizationDirection().Dot(incidentKVector);
    	double incidentCrossPol = currentElement.GetCrossPolarizationDirection().Dot(incidentKVector);


        // theta is elevation: sweeping xy plane
        // phi is azimuth: sweeping xz plane
        // the angle is measured between the incoming radiation and the normal
        // the range is (-pi/2, pi/2), assuming no radiation coming from behind antenna
        // normal incidence is 0

        double tTheta = LMCConst::Pi()/2.; 
        if (incidentNormal>0.) tTheta = atan(incidentCoPol/incidentNormal);
        double thetaFactor = pow(cos(tTheta), 2.3);

        
        double tPhi = LMCConst::Pi()/2.;
        if (incidentNormal>0.) tPhi = atan(incidentCrossPol/incidentNormal);
        double phiFactor = pow(cos(tPhi), 2.3);

        /*
        // dipole theta dependence, normalized to 1.0 for normal incidence
        // SAVED FOR SLOT ANTENNA:
        // theta is the elevation: sweeping xy plane, ranging from 0 to pi
        // the angle is measured between the incoming radiation and the co-pol direction
        // normal incidence is pi/2, parallel to co-pol is 0, antiparallel to co-pol is pi

        double tTheta = LMCConst::Pi()/2.;
        if (incidentCoPol>0.) tTheta = atan(incidentNormal/fabs(incidentCoPol));
        if (incidentCoPol<0.) tTheta = LMCConst::Pi() - atan(incidentNormal/fabs(incidentCoPol));

        // dipole theta dependence, normalized to 1.0 for normal incidence
        double dipoleThetaFactor = 0.0;
        if (fabs(tTheta)>0.) dipoleThetaFactor = cos((LMCConst::Pi()/2.)*cos(tTheta))/sin(tTheta);

    	// dipole donut pinch from HFSS, normalized to 1.0 for normal incidence.
    	double dipolePhiPinchFactor = cos(tPhi);
        */

        //  isotropic placeholder
        /*
        double thetaFactor = 1.0;
        double phiFactor = 1.0;
        */
        

    	return thetaFactor * phiFactor;
    }


} /* namespace locust */

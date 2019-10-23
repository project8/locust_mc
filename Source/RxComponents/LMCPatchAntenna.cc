/*
 * LMCPatchAntenna.cc
 *
 *  Created on: Mar 1, 2018
 *      Author: nbuzinsky
 */

#include "LMCPatchAntenna.hh"


namespace locust
{

    PatchAntenna::PatchAntenna():
        antennaFactorSpline(antennaFactor.begin(), antennaFactor.end(), lowerBoundFrequency, frequencySpacingSpline ),
        gainSpline(gain.begin(), gain.end(), lowerBoundAngle, angularSpacingSpline),
        copolarizationDirection(0,0,0),
        normalDirection(0,0,0),
        centerPosition(0,0,0),
        incidentElectricField(0,0,0),
        incidentMagneticField(0,0,0),
        instantaneousFrequency(0.),
        timeDelay(0.)
    {
    }

    PatchAntenna::~PatchAntenna()
    {
    }
    LMCThreeVector PatchAntenna::GetPosition()
    {
        return centerPosition;
    }

    LMCThreeVector PatchAntenna::GetNormalDirection()
    {
    	return normalDirection;
    }

    void PatchAntenna::SetIncidentElectricField(const LMCThreeVector &incomingElectricField)
    {
        incidentElectricField = incomingElectricField;
    }

    void PatchAntenna::SetIncidentMagneticField(const LMCThreeVector &incomingMagneticField)
    {
        incidentMagneticField = incomingMagneticField;
    }
    void PatchAntenna::SetInstantaneousFrequency(const double &dopplerFrequency)
    {
        instantaneousFrequency = dopplerFrequency;
    }

    double PatchAntenna::GetAnalogTimeDelay()
    {
        return timeDelay;
    }

    double PatchAntenna::GetVoltage()
    {
        return GetCopolarizationFactor() * GetGainFactor() / GetAntennaFactor();
    }

    double PatchAntenna::GetAntennaFactor()
    {
        if(instantaneousFrequency <= 25.1e9 || instantaneousFrequency>=26.1e9)
            return 600;
        return antennaFactorSpline(instantaneousFrequency);
    }

    double PatchAntenna::GetGainFactor()
    {
        LMCThreeVector waveVector = incidentElectricField.Cross(incidentMagneticField);
        waveVector = waveVector.Unit(); //Normalize
        double incidentAngle =  acos(waveVector.Dot(normalDirection));
        if(incidentAngle > LMCConst::Pi() / 2.)
            incidentAngle=LMCConst::Pi() - incidentAngle;

        return  gainSpline(0.) / gainSpline(incidentAngle);
    }

    double PatchAntenna::GetCopolarizationFactor()
    {
      return incidentElectricField.Dot(copolarizationDirection);
    }

    void PatchAntenna::SetCenterPosition(const LMCThreeVector &newPosition)
    {
        centerPosition = newPosition;
    }
    void PatchAntenna::SetPolarizationDirection(const LMCThreeVector &copolDirection)
    {
        copolarizationDirection = copolDirection;
    }

    LMCThreeVector PatchAntenna::GetPolarizationDirection()
     {
         return copolarizationDirection;
     }

    void PatchAntenna::SetNormalDirection(const LMCThreeVector &normDirection)
    {
        normalDirection = normDirection;
    }




} /* namespace locust */

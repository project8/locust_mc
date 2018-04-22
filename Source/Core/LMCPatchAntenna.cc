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
        previousRetardedIndex(-99),
        previousRetardedTime(-99),
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

    int PatchAntenna::GetPreviousRetardedIndex()
    {
        return previousRetardedIndex;
    }

    double PatchAntenna::GetPreviousRetardedTime()
    {
        return previousRetardedTime;
    }

    double PatchAntenna::GetAnalogTimeDelay()
    {
        return timeDelay;
    }

    void PatchAntenna::SetPreviousRetardedIndex(const int& index)
    {
        previousRetardedIndex = index;
    }
    void PatchAntenna::SetPreviousRetardedTime(const double &time)
    {
        previousRetardedTime = time;
    }


    double PatchAntenna::GetVoltage()
    {
        return GetCopolarizationFactor() * GetGainFactor() / GetAntennaFactor();
    }

    double PatchAntenna::GetAntennaFactor()
    {
        return antennaFactorSpline(instantaneousFrequency);
    }

    double PatchAntenna::GetGainFactor()
    {
        LMCThreeVector waveVector = incidentElectricField.Cross(incidentMagneticField);
        waveVector = waveVector.Unit(); //Normalize
        double incidentAngle =  acos(waveVector.Dot(normalDirection));

        return sqrt( gainSpline(incidentAngle)/gainSpline(0.) );
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



} /* namespace locust */

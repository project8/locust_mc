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
        incidentElectricField(0,0,0),
        incidentMagneticField(0,0,0),
        instantaneousFrequency(0.),
        timeDelay(0.)
    {
    }

    PatchAntenna::~PatchAntenna()
    {
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
        double incidentAngle =  acos(waveVector.Dot(GetNormalDirection()));
        if(incidentAngle > LMCConst::Pi() / 2.)
            incidentAngle=LMCConst::Pi() - incidentAngle;

        return  gainSpline(0.) / gainSpline(incidentAngle);
    }

    double PatchAntenna::GetCopolarizationFactor()
    {
      return incidentElectricField.Dot(GetPolarizationDirection());
    }


    void PatchAntenna::RxSayHello()
     {
     	printf("patch says hello\n");
     	getchar();
     }




} /* namespace locust */

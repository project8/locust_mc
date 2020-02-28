/*
 * LMCPlaneWaveTransmitter.hh
 *
 *  Created on: Jan 27, 2020
 *      Author: pslocum
 */

#ifndef LMCPLANEWAVETRANSMITTER_HH_
#define LMCPLANEWAVETRANSMITTER_HH_

#include "LMCTransmitter.hh"
#include "LMCConst.hh"
#include "param.hh"
#include "LMCThreeVector.hh"

namespace locust
{

    /*!
     @class PlaneWaveTransmitter
     @author P. L. Slocum Jan. 27 2020

     @brief Class to generate tranmistter from plane wave

     @details
     Operates in time space

     Configuration name: "plane-wave"

     Available configuration options:

     */
    class PlaneWaveTransmitter : public Transmitter
    {
    public:

        PlaneWaveTransmitter();
        virtual ~PlaneWaveTransmitter();

        bool Configure( const scarab::param_node& aNode );

        virtual double* GetEFieldCoPol(int fieldIndexPoint, double dt);


    private:

        double fAOI; // from json file, in degrees.
        double fAmplitude;
        double fRF_Frequency;  // typically defined by a parameter in json file.
        double fPhaseDelay=0.0; //Delay in the phase that changes for each time sample
        LMCThreeVector fIncidentKVector;  // vector pointing from plane wave to requested point of interest.
        void AddIncidentKVector(LMCThreeVector pointOfInterest);
        void AddPropagationPhaseDelay(LMCThreeVector pointOfInterest);

    };


} /* namespace locust */

#endif /* LMCPLANEWAVETRANSMITTER_HH_ */

/*
 * LMCTransmitter.hh
 *
 *  Created on: Jan. 24, 2020
 *      Author: pslocum
 */

#ifndef LMCTRANSMITTER_HH_
#define LMCTRANSMITTER_HH_

#include "LMCThreeVector.hh"
#include "param.hh"

#include "logger.hh"

#include <vector>

namespace locust
{
 /*!
 @class Transmitter
 @author P. Slocum
 @brief Base class to characterize Transmitter selection
 @details
 Available configuration options:
 No input parameters
 */


    class Transmitter
    {

        public:
            Transmitter();
            virtual ~Transmitter();
            virtual void TxSayHello();

            virtual bool Configure( const scarab::param_node& ){return true;};
            virtual std::vector<double> GetEFieldCoPol(int fieldPointIndex, double dt) {return {0.};};

            virtual std::vector<double> SolveKassFields(LMCThreeVector pointOfInterest, LMCThreeVector coPolDirection, double tReceiverTime, unsigned tTotalElementIndex) {return {0.};};
            virtual std::vector<double> ExtractParticleXP(double TOld, bool Interpolate, double dt, bool EGun) {return {0.};};

            virtual void InitializeFieldPoint(LMCThreeVector fieldPoint);

            virtual bool IsKassiopeia() {return false;};

            /// Initialize the FIR filter and the field estimator
            virtual bool InitializeTransmitter(){return true;};
            LMCThreeVector GetFieldPoint(int index);
            virtual LMCThreeVector GetIncidentKVector(int index);
            double GetPropagationPhaseDelay(int index);
            double GetNPoints();
            double GetMeanofFieldPoints(int axis);
            LMCThreeVector GetMeanofFieldPoints();

	protected:

    	    virtual void AddIncidentKVector(LMCThreeVector fieldPoint);
    	    void AddPropagationPhaseDelay(double phaseDelay);
    	    virtual void AddPropagationPhaseDelay(LMCThreeVector fieldPoint);
    	    void SetFieldPoint(int index,LMCThreeVector fieldPoint);
    	    void SetIncidentKVector(int index,LMCThreeVector incidentKVector);
    	    void SetPropagationPhaseDelay(int index,double phaseDelay);

    private:
    	    std::vector< LMCThreeVector> fFieldPoints;
    	    std::vector< LMCThreeVector> fIncidentKVectors;
    	    std::vector< double> fPropagationPhaseDelays;
    	    void AddFieldPoint(LMCThreeVector fieldPoint);
};


} /* namespace locust */

#endif /* LMCCHANNEL_HH_ */

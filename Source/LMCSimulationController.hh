/*
 * LMCSimulationController.hh
 *
 *  Created on: Feb 5, 2014
 *      Author: nsoblath
 */

#ifndef LMCSIMULATIONCONTROLLER_HH_
#define LMCSIMULATIONCONTROLLER_HH_

#include "LMCEggWriter.hh"
#include "LMCRunLengthCalculator.hh"

namespace locust
{
    class Generator;
    class ParamNode;

    class SimulationController
    {
        public:
            SimulationController();
            virtual ~SimulationController();

            bool Configure( const ParamNode* aNode );

            bool Prepare();

            bool Run() const;

            bool Finalize();

            void SetFirstGenerator( const Generator* firstGen );

        private:
            const Generator* fFirstGenerator;

            RunLengthCalculator fRunLengthCalc;

            EggWriter fEggWriter;
    };

} /* namespace locust */

#endif /* LMCSIMULATIONCONTROLLER_HH_ */

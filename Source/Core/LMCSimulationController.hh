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


namespace scarab
{
    class param_node;
}

namespace locust
{
    class Digitizer;
    class Generator;

    /*!
     @class SimulationController
     @author N. S. Oblath

     @brief Creates and configures the requested generators

     @details

     Configuration name: "simulation"

     Available configuration options:
     - "rng-seed": int -- sets the RNG seed
     - All configuration options in RunLengthCalculator
     - All configuration options in EggWriter

    */
    class SimulationController
    {
        public:
            SimulationController();
            virtual ~SimulationController();

            bool Configure( const scarab::param_node* aNode );

            bool Prepare();

            bool Run();

            bool Finalize();

            void SetRNGSeed(int seed);
            void SetRNGSeed();

            unsigned GetNChannels() const;


            void SetFirstGenerator( Generator* firstGen );

        private:
            const Digitizer* FindDigitizer() const;

            Generator* fFirstGenerator;

            RunLengthCalculator fRunLengthCalc;

            EggWriter fEggWriter;

            const unsigned fNChannels;

    };


    inline unsigned SimulationController::GetNChannels() const
    {
        return fNChannels;
    }



} /* namespace locust */

#endif /* LMCSIMULATIONCONTROLLER_HH_ */

/*
 * LMCDigitizer.hh
 *
 *  Created on: Mar 3, 2014
 *      Author: nsoblath
 */

#ifndef LMCDIGITIZER_HH_
#define LMCDIGITIZER_HH_

#include "LMCGenerator.hh"

#include "thorax.hh"

namespace locust
{

    class Digitizer : public Generator
    {
        public:
            Digitizer( const std::string& aName = "generic-digitizer" );
            virtual ~Digitizer();

            virtual void Configure( const ParamNode* aNode );

            virtual void Accept( GeneratorVisitor* aVisitor ) const;

            const dig_calib_params& DigitizerParams() const;
            dig_calib_params& DigitizerParams();

        protected:
            virtual void Generate( Signal* aSignal ) const;

            struct dig_calib_params fParams;
    };

} /* namespace locust */

#endif /* LMCDIGITIZER_HH_ */

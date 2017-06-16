/*
 * LMCDigitizer.hh
 *
 *  Created on: Mar 3, 2014
 *      Author: nsoblath
 */

#ifndef LMCDIGITIZER_HH_
#define LMCDIGITIZER_HH_

#include "LMCGenerator.hh"

#include "digital.hh"

namespace locust
{

    /*!
     @class Digitizer
     @author N. S. Oblath

     @brief Digitize the data

     @details

     Configuration name: "digitizer"

     Available configuration options:
     - "bit-depth": unsigned -- Number of bits used to digitize the data
     - "data-type-size": unsigned -- Number of bytes used to store the data
     - "v-min": double -- Minimum voltage for the digitizer input
     - "v-range": double -- Range of the digitizer input

    */
    class Digitizer : public Generator
    {
        public:
            Digitizer( const std::string& aName = "digitizer" );
            virtual ~Digitizer();

            bool Configure( const scarab::param_node* aNode );

            void Accept( GeneratorVisitor* aVisitor ) const;

            const scarab::dig_calib_params& DigitizerParams() const;
            scarab::dig_calib_params& DigitizerParams();

            void SetADCValuesSigned( bool aFlag );
            bool GetADCValuesSigned() const;

        protected:
            bool DoGenerate( Signal* aSignal ) const;

            struct scarab::dig_calib_params fParams;

            bool fADCValuesSigned;
    };

    inline void Digitizer::SetADCValuesSigned( bool aFlag )
    {
        fADCValuesSigned = aFlag;
        return;
    }

    inline bool Digitizer::GetADCValuesSigned() const
    {
        return fADCValuesSigned;
    }

} /* namespace locust */

#endif /* LMCDIGITIZER_HH_ */

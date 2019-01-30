
#ifndef LMCPOWERCOMBINER_HH_
#define LMCPOWERCOMBINER_HH_

#include "LMCConst.hh"

namespace locust
{
 /*!
 @class PowerCombiner
 @author P. Slocum
 @brief Class to describe the power combining in the patch array.
 @details
 Available configuration options:
 No input parameters
 */
    class PowerCombiner
    {

        public:
            PowerCombiner();

            virtual ~PowerCombiner();

            double GetLinePhaseCorr(unsigned z_index, double DopplerFrequency);

            double GetVoltageDamping(int njunctions);


        private:

    };


} /* namespace locust */

#endif /* LMCPOWERCOMBINER_HH_ */


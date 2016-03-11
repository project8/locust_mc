/*
 * LMCRunKassiopeia.hh
 *
 *  Created on: Mar 10, 2016
 *      Author: nsoblath
 */

#ifndef LMC_LMCRUNKASSIOPEIA_HH_
#define LMC_LMCRUNKASSIOPEIA_HH_

namespace katrin
{
    class KCommandLineTokenizer;
}

namespace locust
{

    class RunKassiopeia
    {
        public:
            RunKassiopeia( katrin::KCommandLineTokenizer& aCommandLine );
            virtual ~RunKassiopeia();

            int Run();
    };

} /* namespace locust */

#endif /* LMC_LMCRUNKASSIOPEIA_HH_ */

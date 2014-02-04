/*
 * LMCDefaultConfig.hh
 *
 *  Created on: Nov 4, 2013
 *      Author: nsoblath
 */

#ifndef LMCDEFAULTCONFIG_HH_
#define LMCDEFAULTCONFIG_HH_

#include "LMCParam.hh"

namespace locust
{

    class DefaultConfig : public ParamNode
    {
        public:
            DefaultConfig();
            virtual ~DefaultConfig();
    };

} /* namespace locust */
#endif /* LMCDEFAULTCONFIG_HH_ */

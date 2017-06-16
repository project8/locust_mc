/*
 * LMCDefaultConfig.hh
 *
 *  Created on: Nov 4, 2013
 *      Author: nsoblath
 */

#ifndef LMCDEFAULTCONFIG_HH_
#define LMCDEFAULTCONFIG_HH_

#include "param.hh"

namespace locust
{

    class DefaultConfig : public scarab::param_node
    {
        public:
            DefaultConfig();
            virtual ~DefaultConfig();
    };

} /* namespace locust */
#endif /* LMCDEFAULTCONFIG_HH_ */

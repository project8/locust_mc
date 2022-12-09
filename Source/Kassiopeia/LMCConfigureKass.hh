/*
 * LMCConfigureKass.hh
 *
 *  Created on: Oct 20, 2022
 *      Author: pslocum
 */

#ifndef LMCCONFIGUREKASS_HH_
#define LMCCONFIGUREKASS_HH_

#include "param.hh"

namespace locust
{
 /*!
 @class ConfigureKass
 @author P. Slocum
 @brief Class to configure Kassiopeia classes using Scarab parameters.
 @details
 Available configuration options:
 No input parameters
 */
    class ConfigureKass
    {

        public:
            ConfigureKass();
            virtual ~ConfigureKass();
            const scarab::param_node* GetParameters();
            void SetParameters( const scarab::param_node& aNode );


        private:
            const scarab::param_node* fParam;

};


} /* namespace locust */

#endif /* LMCCONFIGURE_HH_ */

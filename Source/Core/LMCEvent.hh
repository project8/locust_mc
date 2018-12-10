/*
 * LMCEvent.hh
 *
 * This class has access to both locust and ROOT libraries.  The syntax is
 * consistent with KTROOTData.hh and the instructions in
 * https://root.cern.ch/root/Using.html .  It is also mentioned in LMCEventLinkDef.hh .
 *  Created on: Dec 5, 2018
 *      Author: pslocum
 */



#ifndef LMCEVENT_HH_
#define LMCEVENT_HH_

#include "TObject.h"


namespace locust
{

    class Event : public TObject 
    {

        public:
            Event();
            virtual ~Event();

            ClassDef(Event,1)  // Root syntax.

    };

}
#endif /* LMCEVENT_HH_ */

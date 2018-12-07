/*
 * LMCEvent.hh
 *
 *  Created on: Dec 5, 2018
 *      Author: pslocum
 */



#ifndef LMCEVENT_HH_
#define LMCEVENT_HH_

#include "TObject.h"


namespace locust
{

    class Event;
    class Event : public TObject 
    {

        public:
            Event();
            virtual ~Event();

   ClassDef(Event,2)  // Root syntax.

};

}
#endif /* LMCEVENT_HH_ */

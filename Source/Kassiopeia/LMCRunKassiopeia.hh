/*
 * LMCRunKassiopeia.hh
 *
 *  Created on: Mar 10, 2016
 *      Author: nsoblath
 */

#ifndef LMC_LMCRUNKASSIOPEIA_HH_
#define LMC_LMCRUNKASSIOPEIA_HH_

#include <vector>
#include <string>

//#include "KElementProcessor.hh"

/*
namespace katrin
{
    class KConditionProcessor;
    class KIncludeProcessor;
    class KLoopProcessor;
    class KPrintProcessor;
    class KTagProcessor;
    class KVariableProcessor;
    class KXMLTokenizer;
    class KCommandLineTokenizer;
//#ifdef Kommon_USE_ROOT
    class KFormulaProcessor;
    class KSaveSettingsProcessor;
//#endif
}
*/


namespace locust
{

    class RunKassiopeia
    {
        public:
            RunKassiopeia();
            virtual ~RunKassiopeia();

            int Run( const std::vector< std::string >& aFiles );
            int Run( const std::string& aFile );

        private:
            /*
            katrin::KXMLTokenizer* fTokenizer;
            katrin::KCommandLineTokenizer* fCommandLineTokenizer;
            katrin::KVariableProcessor* fVariableProcessor;
            katrin::KIncludeProcessor* fIncludeProcessor;
            katrin::KLoopProcessor* fLoopProcessor;
            katrin::KConditionProcessor* fConditionProcessor;
            katrin::KPrintProcessor* fPrintProcessor;
            katrin::KTagProcessor* fTagProcessor;
            katrin::KElementProcessor* fElementProcessor;

//#ifdef Kommon_USE_ROOT
            katrin::KFormulaProcessor* fFormulaProcessor;
//#endif
*/


    };

} /* namespace locust */

#endif /* LMC_LMCRUNKASSIOPEIA_HH_ */


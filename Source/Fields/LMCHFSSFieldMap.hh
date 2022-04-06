/*
 * LMCHFSSFieldMap.hh
 *
 *  Created on: Nov. 17, 2021
 *      Author: pslocum
 */

#ifndef LMCHFSSFIELDMAP_HH_
#define LMCHFSSFIELDMAP_HH_

#include "param.hh"

#include "logger.hh"
#include "LMCField.hh"
#include "LMCKassLocustInterface.hh"

#include <vector>

namespace locust
{
 /*!
 @class HFSSFieldMap
 @author P. Slocum
 @brief Derived class to define HFSSFieldMap fields.
 @details Placeholder class to interpolate mode field maps from HFSS file.
 	 Grid size should be fInterface->fnPixels X fInterface->fnPixels.
 Available configuration options:  Should include filename with field map.
 No input parameters
 */


    class HFSSFieldMap : public Field
    {

        public:
            HFSSFieldMap();
            virtual ~HFSSFieldMap();

            virtual bool Configure( const scarab::param_node& ){return true;};

            void ReadHFSSFile();

            std::vector<double> TE_E(int m, int n, double xKass, double yKass) const;
            std::vector<double> TE_H(int m, int n, double xKass, double yKass) const;
            std::vector<double> TM_E(int m, int n, double xKass, double yKass) const;
            std::vector<double> TM_H(int m, int n, double xKass, double yKass) const;
            double Integrate(int l, int m, int n, bool teMode, bool eField);


        private:
            kl_interface_ptr_t fInterface;


    };


}; /* namespace locust */

#endif /* LMCHFSSFIELDMAP_HH_ */

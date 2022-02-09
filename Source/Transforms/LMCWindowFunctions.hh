/*
 * LMCWindowFunctions.hh
 *
 *  Created on: Jan 10, 2022
 *      Author: A. B. Telles
 */

#ifndef LMCWINDOWFUNCTIONS_HH_
#define LMCWINDOWFUNCTIONS_HH_

#include "param.hh"

namespace locust
{
    /*!
     @class WindowFunctions
     @author A. B. Telles
     @brief Class to create and apply window functions
     @details
     Currently being used for the complex FFT class. Can be expanded to other things such as LPF in the future.
     
     Available configuration options:
     
     */
    class WindowFunctions
    {
        
    public:
        WindowFunctions();
        virtual ~WindowFunctions();
        bool Configure( const scarab::param_node& aNode);
        bool GenerateWindowFunction(int windowsize, int prewindowbins, int signalsize);
        bool SetupWindow(std::string windowname, double windowparam);
        std::vector<double>* GetWindowFunction();
        bool IsWindowGenerated();

    private:
    
        int fWindowFunctionType;
        int fTotalWindowSize;
        int fPreWindowBins;
        int fSignalSize;
        double fTukeyWindowAlpha;
        std::vector<double> fWindowFunction;
        bool IsWindowGeneratedFlag;

        // Member functions:
        void GenerateRectangularWindow();
        void GenerateTukeyWindow();
        void GenerateHanningWindow();
    };
    
} /* namespace locust */

#endif /* LMCWINDOWFUNCTION_HH_ */

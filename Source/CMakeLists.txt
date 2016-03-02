set( LOCUST_MC_HEADER_FILES
    Core/LMCConfigurator.hh
    Core/LMCGlobals.hh
#    Core/LMCConstants.hh
#    Core/LMCDefaultConfig.hh
    Core/LMCDestroyer.hh
 #   Core/LMCEggWriter.hh
    Core/LMCException.hh
    Core/LMCFactory.hh
    Core/LMCGenerator.hh
    Core/LMCGeneratorToolbox.hh
    Core/LMCLogger.hh
    Core/LMCParam.hh
    Core/LMCParser.hh
    Core/LMCRunLengthCalculator.hh
    Core/LMCSignal.hh
 #   Core/LMCSimulationController.hh
    Core/LMCSingleton.hh
    Core/LMCVisitor.hh
    
    Generators/LMCDigitizer.hh
 #   Generators/LMCGaussianNoiseGenerator.hh
    Generators/LMCTestSignalGenerator.hh
    Generators/LMCBasebandTrackGenerator.hh
    Generators/LMCReceiverTransferFunctionGenerator.hh
    Generators/LMCTrappedElectronGenerator.hh
#    Generators/LMCKassiopeiaGenerator.hh
#    Generators/LMCSingleAmpSim.hh
)


set( LOCUST_MC_SOURCE_FILES
    Core/LMCConfigurator.cc
#    Core/LMCDefaultConfig.cc
 #   Core/LMCEggWriter.cc
    Core/LMCException.cc
    Core/LMCGenerator.cc
    Core/LMCGeneratorToolbox.cc
    Core/LMCLogger.cc
    Core/LMCParam.cc
    Core/LMCParser.cc
    Core/LMCRunLengthCalculator.cc
    Core/LMCSignal.cc
 #   Core/LMCSimulationController.cc
    Core/LMCVisitor.cc

    Generators/LMCDigitizer.cc
 #   Generators/LMCGaussianNoiseGenerator.cc
    Generators/LMCTestSignalGenerator.cc
    Generators/LMCBasebandTrackGenerator.cc
    Generators/LMCReceiverTransferFunctionGenerator.cc
    Generators/LMCTrappedElectronGenerator.cc
#    Generators/LMCKassiopeiaGenerator.cc
#    Generators/LMCSingleAmpSim.cc
)


pbuilder_library( LocustMC LOCUST_MC_SOURCE_FILES "" )

pbuilder_install_headers( ${LOCUST_MC_HEADER_FILES} )

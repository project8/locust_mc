set( LOCUST_MC_HEADER_FILES
    Core/LMCConfigurator.hh
    Core/LMCGlobals.hh
#    Core/LMCConstants.hh
#    Core/LMCDefaultConfig.hh
    Core/LMCEggWriter.hh
    Core/LMCException.hh
    Core/LMCGenerator.hh
    Core/LMCGeneratorToolbox.hh
    Core/LMCParam.hh
    Core/LMCParser.hh
    Core/LMCRunLengthCalculator.hh
    Core/LMCSignal.hh
    Core/LMCSimulationController.hh
    Core/LMCVisitor.hh
    
    Generators/LMCDigitizer.hh
    Generators/LMCGaussianNoiseGenerator.hh
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
    Core/LMCEggWriter.cc
    Core/LMCException.cc
    Core/LMCGenerator.cc
    Core/LMCGeneratorToolbox.cc
    Core/LMCParam.cc
    Core/LMCParser.cc
    Core/LMCRunLengthCalculator.cc
    Core/LMCSignal.cc
    Core/LMCSimulationController.cc
    Core/LMCVisitor.cc

    Generators/LMCDigitizer.cc
    Generators/LMCGaussianNoiseGenerator.cc
    Generators/LMCTestSignalGenerator.cc
    Generators/LMCBasebandTrackGenerator.cc
    Generators/LMCReceiverTransferFunctionGenerator.cc
    Generators/LMCTrappedElectronGenerator.cc
#    Generators/LMCKassiopeiaGenerator.cc
#    Generators/LMCSingleAmpSim.cc
)


IF ( locust_mc_BUILD_WITH_KASSIOPEIA )

set( LOCUST_MC_HEADER_FILES
    Core/LMCConfigurator.hh
    Core/LMCGlobals.hh
#    Core/LMCConstants.hh
#    Core/LMCDefaultConfig.hh
    Core/LMCEggWriter.hh
    Core/LMCException.hh
    Core/LMCGenerator.hh
    Core/LMCGeneratorToolbox.hh
    Core/LMCParam.hh
    Core/LMCParser.hh
    Core/LMCRunLengthCalculator.hh
    Core/LMCSignal.hh
    Core/LMCSimulationController.hh
    Core/LMCVisitor.hh
    
    Generators/LMCDigitizer.hh
    Generators/LMCGaussianNoiseGenerator.hh
    Generators/LMCTestSignalGenerator.hh
    Generators/LMCBasebandTrackGenerator.hh
    Generators/LMCReceiverTransferFunctionGenerator.hh
    Generators/LMCTrappedElectronGenerator.hh
    Generators/LMCKassiopeiaGenerator.hh
#    Generators/LMCSingleAmpSim.hh
)


set( LOCUST_MC_SOURCE_FILES
    Core/LMCConfigurator.cc
#    Core/LMCDefaultConfig.cc
    Core/LMCEggWriter.cc
    Core/LMCException.cc
    Core/LMCGenerator.cc
    Core/LMCGeneratorToolbox.cc
    Core/LMCParam.cc
    Core/LMCParser.cc
    Core/LMCRunLengthCalculator.cc
    Core/LMCSignal.cc
    Core/LMCSimulationController.cc
    Core/LMCVisitor.cc

    Generators/LMCDigitizer.cc
    Generators/LMCGaussianNoiseGenerator.cc
    Generators/LMCTestSignalGenerator.cc
    Generators/LMCBasebandTrackGenerator.cc
    Generators/LMCReceiverTransferFunctionGenerator.cc
    Generators/LMCTrappedElectronGenerator.cc
    Generators/LMCKassiopeiaGenerator.cc
#    Generators/LMCSingleAmpSim.cc
)



ENDIF ( locust_mc_BUILD_WITH_KASSIOPEIA )


pbuilder_library( LocustMC LOCUST_MC_SOURCE_FILES "" )
SET_TARGET_PROPERTIES(LocustMC PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
#target_link_libraries( LocustMC ${EXTERNAL_LIBRARIES} )
#target_link_libraries( LocustMC ${KommonCore} )
#target_link_libraries( LocustMC ${KommonRoot} )
#target_link_libraries( LocustMC ${KommonGsl} )
#target_link_libraries( LocustMC ${KommonVtk} )
#target_link_libraries( LocustMC ${KGeoBagCore} )
#target_link_libraries( LocustMC ${KEMCore} )
pbuilder_install_headers( ${LOCUST_MC_HEADER_FILES} )



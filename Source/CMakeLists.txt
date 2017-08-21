
if (locust_mc_BUILD_WITH_KASSIOPEIA)

    # Make up for variables not created by PackageBuilder, but required by kasper_module_paths
    set_path(DATA_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/data" "Install directory for data files")
    set_path(SCRATCH_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/scratch" "Directory for temporary files")
    set_path(OUTPUT_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/output" "Directory for output files")
    set_path(CACHE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/cache" "Directory for cache files")
    set_path(LOG_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/log" "Directory for log files")
    
    kasper_module_paths( locust_mc )
    
    foreach( FLAG ${Kassiopeia_CFLAGS} )
        message( STATUS "%%% adding {FLAG}" )
        add_definitions( -D${FLAG} )
    endforeach( FLAG )
        
    # include directories
    #kasper_internal_include_directories(
    #)

endif (locust_mc_BUILD_WITH_KASSIOPEIA)


set( LOCUST_MC_HEADER_FILES
    #Core/LMCConfigurator.hh
    Core/LMCGlobalsDeclaration.hh
    Core/LMCGlobalsDefinition.hh
#    Core/LMCConstants.hh
#    Core/LMCDefaultConfig.hh
    Core/LMCEggWriter.hh
    Core/LMCException.hh
    Core/LMCGenerator.hh
    Core/LMCGeneratorToolbox.hh
    Core/LMCHFSSReader.hh
    #Core/LMCParam.hh
    #Core/LMCParser.hh
    Core/LMCParticle.hh
    Core/LMCRunLengthCalculator.hh
    Core/LMCSignal.hh
    Core/LMCSimulationController.hh
    Core/LMCVisitor.hh
    
    Generators/LMCDigitizer.hh
    Generators/LMCGaussianNoiseGenerator.hh
#    Generators/LMCTestSignalGenerator.hh
#    Generators/LMCBasebandTrackGenerator.hh
#    Generators/LMCTrappedElectronGenerator.hh
#    Generators/LMCSingleAmpSim.hh
)


set( LOCUST_MC_SOURCE_FILES
    #Core/LMCConfigurator.cc
#    Core/LMCDefaultConfig.cc
    Core/LMCEggWriter.cc
    Core/LMCException.cc
    Core/LMCGenerator.cc
    Core/LMCGeneratorToolbox.cc
    Core/LMCHFSSReader.cc
    #Core/LMCParam.cc
    #Core/LMCParser.cc
    Core/LMCParticle.cc
    Core/LMCRunLengthCalculator.cc
    Core/LMCSignal.cc
    Core/LMCSimulationController.cc
    Core/LMCVisitor.cc

    Generators/LMCDigitizer.cc
    Generators/LMCGaussianNoiseGenerator.cc
#    Generators/LMCTestSignalGenerator.cc
#    Generators/LMCBasebandTrackGenerator.cc
#    Generators/LMCTrappedElectronGenerator.cc
#    Generators/LMCSingleAmpSim.cc
)

if (locust_mc_BUILD_WITH_KASSIOPEIA)
    set( LOCUST_MC_HEADER_FILES ${LOCUST_MC_HEADER_FILES}
        Generators/LMCFreeFieldSignalGenerator.hh
        Generators/LMCDecimateSignalGenerator.hh
        Generators/LMCLowPassFilterFFTGenerator.hh
        Generators/LMCButterworthLPFGenerator.hh


        Kassiopeia/LMCCyclotronRadiationExtractor.hh
        Kassiopeia/LMCCyclotronRadiationExtractorBuilder.hh
        Kassiopeia/LMCEventHold.hh
        Kassiopeia/LMCEventHoldBuilder.hh
        Kassiopeia/LMCRootBuilder.hh
        Kassiopeia/LMCRunKassiopeia.hh
    )

    set( LOCUST_MC_SOURCE_FILES ${LOCUST_MC_SOURCE_FILES}
        Generators/LMCFreeFieldSignalGenerator.cc
        Generators/LMCDecimateSignalGenerator.cc
        Generators/LMCLowPassFilterFFTGenerator.cc
        Generators/LMCButterworthLPFGenerator.cc


        Kassiopeia/LMCCyclotronRadiationExtractor.cc
        Kassiopeia/LMCCyclotronRadiationExtractorBuilder.cc
        Kassiopeia/LMCEventHold.cc
        Kassiopeia/LMCEventHoldBuilder.cc
        Kassiopeia/LMCRootBuilder.cc
        Kassiopeia/LMCRunKassiopeia.cc
    )
endif (locust_mc_BUILD_WITH_KASSIOPEIA)


pbuilder_library( LocustMC LOCUST_MC_SOURCE_FILES "" )
pbuilder_install_headers( ${LOCUST_MC_HEADER_FILES} )



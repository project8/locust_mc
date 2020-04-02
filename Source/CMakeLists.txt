
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
    Core/LMCConst.hh
    Core/LMCEggWriter.hh
    Core/LMCException.hh
    Core/LMCGenerator.hh
    Core/LMCGeneratorToolbox.hh
    Core/LMCRunLengthCalculator.hh
    Core/LMCSignal.hh
    Core/LMCSimulationController.hh
    Core/LMCThreeVector.hh
    Core/LMCVisitor.hh
	Core/LMCHFSSResponseFileHandler.hh
	Core/LMCFIRFileHandler.hh
	Core/LMCTFFileHandler.hh	
    Core/LMCTriangle.cc

    Distributions/LMCBaseDistribution.hh
    Distributions/LMCDiracDistribution.hh
    Distributions/LMCExponentialDistribution.hh
    Distributions/LMCGaussianDistribution.hh
    Distributions/LMCLorentzianDistribution.hh
    Distributions/LMCRuddDistribution.hh
    Distributions/LMCUniformDistribution.hh
    Distributions/LMCKrComplexLineDistribution.hh
    Distributions/LMCDistributionInterface.hh
	    
    Generators/LMCDigitizer.hh
    Generators/LMCGaussianNoiseGenerator.hh
    Generators/LMCTestSignalGenerator.hh
    Generators/LMCDecimateSignalGenerator.hh
    Generators/LMCLowPassFilterFFTGenerator.hh
    Generators/LMCHighPassFilterFFTGenerator.hh
    Generators/LMCLocalOscillatorGenerator.hh
    Generators/LMCButterworthLPFGenerator.hh

    Transmitters/LMCFieldBuffer.hh
   	Transmitters/LMCTransmitter.hh
    Transmitters/LMCAntennaSignalTransmitter.hh
    Transmitters/LMCPlaneWaveTransmitter.hh
    Transmitters/LMCTransmitterHardware.hh
    Transmitters/LMCDipoleAntenna.hh
    Transmitters/LMCTurnstileAntenna.hh
    
    Transforms/LMCHilbertTransform.hh
    Transforms/LMCComplexFFT.hh    

    RxComponents/LMCPowerCombiner.hh
    RxComponents/LMCVoltageDivider.hh
    RxComponents/LMCSlottedWaveguide.hh
    RxComponents/LMCSinglePatch.hh
    RxComponents/LMCCorporateFeed.hh
    RxComponents/LMCsMatrix.hh
    RxComponents/LMCUnitCell.hh
    RxComponents/LMCSeriesFeed.hh
    RxComponents/LMCReceiver.hh
    RxComponents/LMCPatchAntenna.hh
    RxComponents/LMCSlotAntenna.hh
    RxComponents/LMCChannel.hh
    
)


set( LOCUST_MC_SOURCE_FILES
    Core/LMCEggWriter.cc
    Core/LMCException.cc
    Core/LMCGenerator.cc
    Core/LMCGeneratorToolbox.cc
    Core/LMCRunLengthCalculator.cc
    Core/LMCSignal.cc
    Core/LMCSimulationController.cc
    Core/LMCThreeVector.cc
    Core/LMCVisitor.cc
	Core/LMCHFSSResponseFileHandler.cc
    Core/LMCFIRFileHandler.cc
    Core/LMCTFFileHandler.cc
    Core/LMCTriangle.cc

    Distributions/LMCBaseDistribution.cc
    Distributions/LMCDiracDistribution.cc
    Distributions/LMCExponentialDistribution.cc
    Distributions/LMCGaussianDistribution.cc
    Distributions/LMCLorentzianDistribution.cc
    Distributions/LMCRuddDistribution.cc
    Distributions/LMCUniformDistribution.cc
    Distributions/LMCKrComplexLineDistribution.cc
    Distributions/LMCDistributionInterface.cc
    	
    Generators/LMCDigitizer.cc
    Generators/LMCGaussianNoiseGenerator.cc
    Generators/LMCTestSignalGenerator.cc
    Generators/LMCDecimateSignalGenerator.cc
    Generators/LMCLowPassFilterFFTGenerator.cc
    Generators/LMCHighPassFilterFFTGenerator.cc
    Generators/LMCLocalOscillatorGenerator.cc
    Generators/LMCButterworthLPFGenerator.cc
    

	Transmitters/LMCTransmitter.cc
   	Transmitters/LMCAntennaSignalTransmitter.cc
    Transmitters/LMCPlaneWaveTransmitter.cc
   	Transmitters/LMCFieldBuffer.cc
   	Transmitters/LMCTransmitterHardware.cc
   	Transmitters/LMCDipoleAntenna.cc
   	Transmitters/LMCTurnstileAntenna.cc
    
    Transforms/LMCHilbertTransform.cc
    Transforms/LMCComplexFFT.cc    
    
    RxComponents/LMCPowerCombiner.cc
    RxComponents/LMCVoltageDivider.cc
    RxComponents/LMCSlottedWaveguide.cc
    RxComponents/LMCSinglePatch.cc
    RxComponents/LMCCorporateFeed.cc
    RxComponents/LMCsMatrix.cc
    RxComponents/LMCUnitCell.cc
    RxComponents/LMCSeriesFeed.cc
    RxComponents/LMCReceiver.cc
    RxComponents/LMCPatchAntenna.cc
    RxComponents/LMCSlotAntenna.cc
    RxComponents/LMCChannel.cc

)

if (locust_mc_BUILD_WITH_KASSIOPEIA)
    set( LOCUST_MC_HEADER_FILES ${LOCUST_MC_HEADER_FILES}

        Generators/LMCFreeFieldSignalGenerator.hh
        Generators/LMCKassSignalGenerator.hh
	Generators/LMCTransmitterInterfaceGenerator.hh
	Generators/LMCFieldParameterGenerator.hh
	Generators/LMCArraySignalGenerator.hh

        Transmitters/LMCKassTransmitter.hh

        Core/LMCGlobalsDeclaration.hh
        Core/LMCGlobalsDefinition.hh
        Core/LMCParticle.hh
        Core/LMCLienardWiechert.hh
        Core/LMCHFSSReader.hh

        Kassiopeia/LMCCyclotronRadiationExtractor.hh
        Kassiopeia/LMCCyclotronRadiationExtractorBuilder.hh
        Kassiopeia/LMCFieldCalculator.hh
        Kassiopeia/LMCFieldCalculatorBuilder.hh
        Kassiopeia/LMCEventHold.hh
        Kassiopeia/LMCEventHoldBuilder.hh
        Kassiopeia/LMCRootBuilder.hh
        Kassiopeia/LMCRunKassiopeia.hh
    )

    set( LOCUST_MC_SOURCE_FILES ${LOCUST_MC_SOURCE_FILES}
        Generators/LMCFreeFieldSignalGenerator.cc
        Generators/LMCKassSignalGenerator.cc
	Generators/LMCTransmitterInterfaceGenerator.cc
	Generators/LMCFieldParameterGenerator.cc
	Generators/LMCArraySignalGenerator.cc        
        
        Transmitters/LMCKassTransmitter.cc
        
        Core/LMCHFSSReader.cc
        Core/LMCParticle.cc
        Core/LMCLienardWiechert.cc

        Kassiopeia/LMCFieldCalculator.cc
        Kassiopeia/LMCFieldCalculatorBuilder.cc
        Kassiopeia/LMCEventHold.cc
        Kassiopeia/LMCEventHoldBuilder.cc
        Kassiopeia/LMCRootBuilder.cc
        Kassiopeia/LMCRunKassiopeia.cc
    )
endif (locust_mc_BUILD_WITH_KASSIOPEIA)



#  Check root build.
if (locust_mc_BUILD_WITH_ROOT)

    set( LOCUST_MC_HEADER_FILES ${LOCUST_MC_HEADER_FILES}
        Core/LMCEvent.hh
        Core/LMCTrack.hh
        Generators/LMCFakeTrackSignalGenerator.hh
       )

    set( LOCUST_MC_SOURCE_FILES ${LOCUST_MC_SOURCE_FILES}
        Core/LMCEvent.cc
        Core/LMCTrack.cc
        Generators/LMCFakeTrackSignalGenerator.cc
       )


    if (CMAKE_COMPILER_IS_GNUCXX)
        option( GCC_FORCE_LINKING "Fix linker errors with some GCC versions by adding the --no-as-needed flag" ON )
        if( GCC_FORCE_LINKING )
            set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed" )
        endif( )
    endif()

    pbuilder_add_ext_libraries (${ROOT_LIBRARIES})
    message(STATUS "${ROOT_LIBRARIES}")
    
    set (DICT_OUTFILES
        ${CMAKE_CURRENT_BINARY_DIR}/LMCDict.cxx
    )
      
    set (DICT_PCMFILES 
        ${CMAKE_CURRENT_BINARY_DIR}/LMCDict_rdict.pcm
    )

    # root_generate_dictionary is not finding everything.  use add_custom_command instead:
    # run rootcint statements on command line from cbuild directory instead.
    #    ROOT_GENERATE_DICTIONARY(EventDict ${EVENT_DICT_HEADERFILES} LINKDEF ${EVENT_LINKDEF_HEADERFILE})   
    get_directory_property(incdirs INCLUDE_DIRECTORIES)
    set(includedirs)
    foreach( d ${incdirs})
       set(includedirs ${includedirs} -I${d})
    endforeach()
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/LMCDict.cxx 
        COMMAND rootcling -f LMCDict.cxx -c -inlineInputHeader ${includedirs}
            ${PROJECT_SOURCE_DIR}/Source/Core/LMCEvent.hh 
            ${PROJECT_SOURCE_DIR}/Source/Core/LMCTrack.hh 
            ${PROJECT_SOURCE_DIR}/Source/Core/LinkDef/LMCLinkDef.hh 
        VERBATIM
    )
    
    pbuilder_install_files (${LIB_INSTALL_DIR} ${DICT_PCMFILES})
        
endif(locust_mc_BUILD_WITH_ROOT)


#pbuilder_library( LocustMC LOCUST_MC_SOURCE_FILES EVENT_DICT_OUTFILE)
add_library(LocustMC ${LOCUST_MC_SOURCE_FILES} ${DICT_OUTFILES})
target_link_libraries (LocustMC ${EXTERNAL_LIBRARIES})

pbuilder_install_headers( ${LOCUST_MC_HEADER_FILES} )
pbuilder_install_libraries( LocustMC )

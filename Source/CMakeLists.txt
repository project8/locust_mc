
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
    Core/LMCFIRHandler.hh
    
    Generators/LMCDigitizer.hh
    Generators/LMCGaussianNoiseGenerator.hh
    Generators/LMCTestSignalGenerator.hh
    Generators/LMCDecimateSignalGenerator.hh
    Generators/LMCLowPassFilterFFTGenerator.hh
    Generators/LMCHighPassFilterFFTGenerator.hh
    Generators/LMCLocalOscillatorGenerator.hh
    Generators/LMCButterworthLPFGenerator.hh
    Generators/LMCDipoleSignalGenerator.cc
    Generators/LMCPlaneWaveSignalGenerator.hh

    Transmitters/LMCFieldBuffer.hh
    Transmitters/LMCAntennaSignalTransmitter.hh
    Transforms/LMCHilbertTransform.hh
    
    RxComponents/LMCPowerCombiner.hh
    RxComponents/LMCReceiver.hh
    RxComponents/LMCPatchAntenna.hh
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
    Core/LMCFIRHandler.cc

    Generators/LMCDigitizer.cc
    Generators/LMCGaussianNoiseGenerator.cc
    Generators/LMCTestSignalGenerator.cc
    Generators/LMCDecimateSignalGenerator.cc
    Generators/LMCLowPassFilterFFTGenerator.cc
    Generators/LMCHighPassFilterFFTGenerator.cc
    Generators/LMCLocalOscillatorGenerator.cc
    Generators/LMCButterworthLPFGenerator.cc
    Generators/LMCDipoleSignalGenerator.cc
    Generators/LMCPlaneWaveSignalGenerator.cc

    Transmitters/LMCFieldBuffer.cc
    Transmitters/LMCAntennaSignalTransmitter.cc
    Transforms/LMCHilbertTransform.cc
    
    RxComponents/LMCPowerCombiner.cc
    RxComponents/LMCReceiver.cc
    RxComponents/LMCPatchAntenna.cc
    RxComponents/LMCChannel.cc

)

if (locust_mc_BUILD_WITH_KASSIOPEIA)
    # Here we only need to add the Kassiopeia-dependent header/source files to the header/source-file lists
    
    set( LOCUST_MC_HEADER_FILES ${LOCUST_MC_HEADER_FILES}

        #Generators/LMCFreeFieldSignalGenerator.hh
        #Generators/LMCPatchSignalGenerator.hh
        Generators/LMCKassSignalGenerator.hh

        Core/LMCGlobalsDeclaration.hh
        Core/LMCGlobalsDefinition.hh
        Core/LMCParticle.hh
        Core/LMCHFSSReader.hh

        Kassiopeia/LMCCyclotronRadiationExtractor.hh
        Kassiopeia/LMCCyclotronRadiationExtractorBuilder.hh
        Kassiopeia/LMCFieldCalculator.hh
        Kassiopeia/LMCEventHold.hh
        Kassiopeia/LMCEventHoldBuilder.hh
        Kassiopeia/LMCKassLocustInterface.hh
        Kassiopeia/LMCParticle.hh
        Kassiopeia/LMCRunKassiopeia.hh
        Kassiopeia/LMCRunPause.hh
        Kassiopeia/LMCRunPauseBuilder.hh
    )

    set( LOCUST_MC_SOURCE_FILES ${LOCUST_MC_SOURCE_FILES}
        #Generators/LMCFreeFieldSignalGenerator.cc
        #Generators/LMCPatchSignalGenerator.cc        
        Generators/LMCKassSignalGenerator.cc

        # Why is this here?
        Core/LMCHFSSReader.cc

        Kassiopeia/LMCCyclotronRadiationExtractor.cc
        Kassiopeia/LMCCyclotronRadiationExtractorBuilder.cc
        Kassiopeia/LMCFieldCalculator.cc
        Kassiopeia/LMCEventHold.cc
        Kassiopeia/LMCEventHoldBuilder.cc
        Kassiopeia/LMCKassLocustInterface.cc
        Kassiopeia/LMCParticle.cc
        Kassiopeia/LMCRunKassiopeia.cc
        Kassiopeia/LMCRunPause.cc
        Kassiopeia/LMCRunPauseBuilder.cc
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

    set (LOCUST_MC_SOURCE_FILES
        ${LOCUST_MC_SOURCE_FILES}
        ${DICT_OUTFILES}
    )
    
    pbuilder_install_files (${LIB_INSTALL_DIR} ${DICT_PCMFILES})
        
endif(locust_mc_BUILD_WITH_ROOT)


pbuilder_library( LocustMC LOCUST_MC_SOURCE_FILES "" )

pbuilder_install_headers( ${LOCUST_MC_HEADER_FILES} )

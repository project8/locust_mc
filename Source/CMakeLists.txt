
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
#    Core/LMCHFSSReader.hh
#    Core/LMCParticle.hh
    Core/LMCRunLengthCalculator.hh
    Core/LMCSignal.hh
    Core/LMCSimulationController.hh
    Core/LMCThreeVector.hh
    Core/LMCVisitor.hh
    
    Generators/LMCDigitizer.hh
    Generators/LMCGaussianNoiseGenerator.hh
    Generators/LMCFakeTrackSignalGenerator.hh    
    Generators/LMCTestSignalGenerator.hh
    Generators/LMCDecimateSignalGenerator.hh
    Generators/LMCLowPassFilterFFTGenerator.hh
    Generators/LMCLocalOscillatorGenerator.hh
    Generators/LMCButterworthLPFGenerator.hh

#    Generators/LMCBasebandTrackGenerator.hh
#    Generators/LMCTrappedElectronGenerator.hh
)


set( LOCUST_MC_SOURCE_FILES
    Core/LMCEggWriter.cc
    Core/LMCException.cc
    Core/LMCGenerator.cc
    Core/LMCGeneratorToolbox.cc
#    Core/LMCHFSSReader.cc
#    Core/LMCParticle.cc
    Core/LMCRunLengthCalculator.cc
    Core/LMCSignal.cc
    Core/LMCSimulationController.cc
    Core/LMCThreeVector.cc
    Core/LMCVisitor.cc

    Generators/LMCDigitizer.cc
    Generators/LMCGaussianNoiseGenerator.cc
    Generators/LMCFakeTrackSignalGenerator.cc    
    Generators/LMCTestSignalGenerator.cc
    Generators/LMCDecimateSignalGenerator.cc
    Generators/LMCLowPassFilterFFTGenerator.cc
    Generators/LMCLocalOscillatorGenerator.cc
    Generators/LMCButterworthLPFGenerator.cc

#    Generators/LMCBasebandTrackGenerator.cc
#    Generators/LMCTrappedElectronGenerator.cc
)

if (locust_mc_BUILD_WITH_ROOT)
    set( LOCUST_MC_HEADER_FILES ${LOCUST_MC_HEADER_FILES}
        Core/LMCEvent.hh
        Core/LMCTrack.hh)
    set( LOCUST_MC_SOURCE_FILES ${LOCUST_MC_SOURCE_FILES}
        Core/LMCEvent.cc
        Core/LMCTrack.cc)
endif (locust_mc_BUILD_WITH_ROOT)

if (locust_mc_BUILD_WITH_KASSIOPEIA)
    set( LOCUST_MC_HEADER_FILES ${LOCUST_MC_HEADER_FILES}

        Generators/LMCFreeFieldSignalGenerator.hh
        Generators/LMCPatchSignalGenerator.hh
        Generators/LMCKassSignalGenerator.hh

        Core/LMCChannel.hh
        Core/LMCGlobalsDeclaration.hh
        Core/LMCGlobalsDefinition.hh
        Core/LMCPatchAntenna.hh
        Core/LMCParticle.hh
        Core/LMCReceiver.hh
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
        Generators/LMCPatchSignalGenerator.cc        
        Generators/LMCKassSignalGenerator.cc

        Core/LMCChannel.cc
        Core/LMCHFSSReader.cc
        Core/LMCParticle.cc
        Core/LMCPatchAntenna.cc
        Core/LMCReceiver.cc

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

  if (CMAKE_COMPILER_IS_GNUCXX)
      option( GCC_FORCE_LINKING "Fix linker errors with some GCC versions by adding the --no-as-needed flag" ON )
    if( GCC_FORCE_LINKING )
        set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed" )
    endif( )
  endif()

    pbuilder_add_ext_libraries (${ROOT_LIBRARIES})
    message(STATUS "${ROOT_LIBRARIES}")
    
    set (EVENT_DICT_HEADERFILES LMCEvent.hh LMCTrack.hh)
    set (EVENT_LINKDEF_HEADERFILE LinkDef/LMCEventLinkDef.hh LMCTrackLinkDef.hh)
# The next two lines are what should be happening, but without root_generate_dictionary they are not working yet.
#    set (EVENT_DICT_OUTFILE ${CMAKE_CURRENT_BINARY_DIR}/EventDict.cxx)
#    set (EVENT_DICT_PCMFILE ${CMAKE_CURRENT_BINARY_DIR}/EventDict_rdict.pcm)

# starting band-aid lines instead:
    set (EVENT_DICT_OUTFILE 
      ${PROJECT_SOURCE_DIR}/Source/EventDict.cxx
      ${PROJECT_SOURCE_DIR}/Source/TrackDict.cxx     
      )
    set (EVENT_DICT_PCMFILE 
      ${PROJECT_SOURCE_DIR}/Source/EventDict_rdict.pcm
      ${PROJECT_SOURCE_DIR}/Source/TrackDict_rdict.pcm      
      )
# end band-aid.

# root_generate_dictionary is not working due to possible path conflict with Scarab.
# run rootcint statements on command line from cbuild directory instead.
    ROOT_GENERATE_DICTIONARY(EventDict ${EVENT_DICT_HEADERFILES} LINKDEF ${EVENT_LINKDEF_HEADERFILE})   
#    rootcint -f EventDict.cxx -c ../Source/Core/LMCEvent.hh ../Source/Core/LinkDef/LMCEventLinkDef.hh 
#    rootcint -f TrackDict.cxx -c ../Source/Core/LMCTrack.hh ../Source/Core/LinkDef/LMCTrackLinkDef.hh 

    pbuilder_install_files (${LIB_INSTALL_DIR}  ${EVENT_DICT_PCMFILE})
        
endif(locust_mc_BUILD_WITH_ROOT)

pbuilder_library( LocustMC LOCUST_MC_SOURCE_FILES EVENT_DICT_OUTFILE )
pbuilder_install_headers( ${LOCUST_MC_HEADER_FILES} )



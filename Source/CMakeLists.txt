
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

    Distributions/LMCBaseDistribution.hh
    Distributions/LMCDiracDistribution.hh
    Distributions/LMCExponentialDistribution.hh
    Distributions/LMCGaussianDistribution.hh
    Distributions/LMCLorentzianDistribution.hh
    Distributions/LMCRuddDistribution.hh
    Distributions/LMCUniformCircleDistribution.hh
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
    Transmitters/LMCTransmitterInterface.hh
    Transmitters/LMCTransmitter.hh
    Transmitters/LMCAntennaSignalTransmitter.hh
    Transmitters/LMCPlaneWaveTransmitter.hh
    Transmitters/LMCTransmitterHardware.hh
    Transmitters/LMCDipoleAntenna.hh
    Transmitters/LMCTurnstileAntenna.hh
    
    Transforms/LMCHilbertTransform.hh
    Transforms/LMCComplexFFT.hh 
    Transforms/LMCWindowFunctions.hh   

    RxComponents/LMCPowerCombinerInterface.hh
    RxComponents/LMCPowerCombiner.hh
    RxComponents/LMCVoltageDivider.hh
    RxComponents/LMCSlottedWaveguide.hh
    RxComponents/LMCCavityModes.hh
    RxComponents/LMCSinglePatch.hh
    RxComponents/LMCCorporateFeed.hh
    RxComponents/LMCsMatrix.hh
    RxComponents/LMCUnitCell.hh
    RxComponents/LMCSeriesFeed.hh
    RxComponents/LMCReceiver.hh
    RxComponents/LMCPatchAntenna.hh
    RxComponents/LMCSlotAntenna.hh
    RxComponents/LMCChannel.hh
    RxComponents/LMCAntennaElementPositionerInterface.hh
    RxComponents/LMCAntennaElementPositioner.hh
    RxComponents/LMCSinglePatchPositioner.hh
    RxComponents/LMCPlanarArrayPositioner.hh
    RxComponents/LMCAnalyticResponseFunction.hh
    RxComponents/LMCEquivalentCircuit.hh
    RxComponents/LMCDampedHarmonicOscillator.hh
    RxComponents/LMCWaveguide.hh

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

    Distributions/LMCBaseDistribution.cc
    Distributions/LMCDiracDistribution.cc
    Distributions/LMCExponentialDistribution.cc
    Distributions/LMCGaussianDistribution.cc
    Distributions/LMCLorentzianDistribution.cc
    Distributions/LMCRuddDistribution.cc
    Distributions/LMCUniformCircleDistribution.cc
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
    Transforms/LMCWindowFunctions.cc    
    
    RxComponents/LMCPowerCombiner.cc
    RxComponents/LMCVoltageDivider.cc
    RxComponents/LMCSlottedWaveguide.cc
    RxComponents/LMCCavityModes.cc
    RxComponents/LMCSinglePatch.cc
    RxComponents/LMCCorporateFeed.cc
    RxComponents/LMCsMatrix.cc
    RxComponents/LMCUnitCell.cc
    RxComponents/LMCSeriesFeed.cc
    RxComponents/LMCReceiver.cc
    RxComponents/LMCPatchAntenna.cc
    RxComponents/LMCSlotAntenna.cc
    RxComponents/LMCChannel.cc    
    RxComponents/LMCAntennaElementPositioner.cc
    RxComponents/LMCSinglePatchPositioner.cc
    RxComponents/LMCPlanarArrayPositioner.cc
    RxComponents/LMCAnalyticResponseFunction.cc
    RxComponents/LMCEquivalentCircuit.cc
    RxComponents/LMCDampedHarmonicOscillator.cc
    RxComponents/LMCWaveguide.cc


)

if (locust_mc_BUILD_WITH_KASSIOPEIA)
    # Here we only need to add the Kassiopeia-dependent header/source files to the header/source-file lists
    
    set( LOCUST_MC_HEADER_FILES ${LOCUST_MC_HEADER_FILES}

        Generators/LMCFreeFieldSignalGenerator.hh
        Generators/LMCKassSignalGenerator.hh
        Generators/LMCArraySignalGenerator.hh
        Generators/LMCCavitySignalGenerator.hh


        Transmitters/LMCKassTransmitter.hh
        Transmitters/LMCKassCurrentTransmitter.hh
        
        RxComponents/LMCWaveguideModes.hh        
        
        Core/LMCHFSSReader.hh

        Kassiopeia/LMCCyclotronRadiationExtractor.hh
        Kassiopeia/LMCCyclotronRadiationExtractorBuilder.hh
        Kassiopeia/LMCLienardWiechert.hh
        Kassiopeia/LMCParticle.hh
        Kassiopeia/LMCConfigureKass.hh
        Kassiopeia/LMCFieldCalculator.hh
        Kassiopeia/LMCEventHold.hh
        Kassiopeia/LMCEventHoldBuilder.hh
        Kassiopeia/LMCKassLocustInterface.hh
        Kassiopeia/LMCRunKassiopeia.hh
        Kassiopeia/LMCRunPause.hh
        Kassiopeia/LMCRunPauseBuilder.hh
        
        Fields/LMCField.hh
        Fields/LMCCylindricalCavity.hh
        Fields/LMCRectangularWaveguide.hh
        Fields/LMCHFSSFieldMap.hh
    )

    set( LOCUST_MC_SOURCE_FILES ${LOCUST_MC_SOURCE_FILES}
        Generators/LMCFreeFieldSignalGenerator.cc
        Generators/LMCKassSignalGenerator.cc
        Generators/LMCArraySignalGenerator.cc   
        Generators/LMCCavitySignalGenerator.cc   
        
        Transmitters/LMCKassTransmitter.cc
        Transmitters/LMCKassCurrentTransmitter.cc
        
        RxComponents/LMCWaveguideModes.cc        
        
        # Why is this here?
        Core/LMCHFSSReader.cc

        Kassiopeia/LMCCyclotronRadiationExtractor.cc
        Kassiopeia/LMCCyclotronRadiationExtractorBuilder.cc
        Kassiopeia/LMCParticle.cc
        Kassiopeia/LMCConfigureKass.cc
        Kassiopeia/LMCLienardWiechert.cc
        Kassiopeia/LMCFieldCalculator.cc
        Kassiopeia/LMCEventHold.cc
        Kassiopeia/LMCEventHoldBuilder.cc
        Kassiopeia/LMCKassLocustInterface.cc
        Kassiopeia/LMCRunKassiopeia.cc
        Kassiopeia/LMCRunPause.cc
        Kassiopeia/LMCRunPauseBuilder.cc
        
        Fields/LMCField.cc
        Fields/LMCCylindricalCavity.cc
        Fields/LMCRectangularWaveguide.cc
        Fields/LMCHFSSFieldMap.cc
    )
endif (locust_mc_BUILD_WITH_KASSIOPEIA)



#  Check root build.
if (locust_mc_BUILD_WITH_ROOT)

    set( LOCUST_MC_HEADER_FILES ${LOCUST_MC_HEADER_FILES}
        Core/LMCEvent.hh
        Core/LMCTrack.hh
        Core/LMCRunParameters.hh
        Generators/LMCFakeTrackSignalGenerator.hh    
        Generators/LMCFakeFreeSpaceSignalGenerator.hh
        IO/LMCFileWriter.hh
    	IO/LMCRootTreeWriter.hh
    	IO/LMCRootHistoWriter.hh
    	IO/LMCRootGraphWriter.hh
       )

    set( LOCUST_MC_SOURCE_FILES ${LOCUST_MC_SOURCE_FILES}
        Core/LMCEvent.cc
        Core/LMCTrack.cc
        Core/LMCRunParameters.cc
        Generators/LMCFakeTrackSignalGenerator.cc   
        Generators/LMCFakeFreeSpaceSignalGenerator.cc
        IO/LMCFileWriter.cc
    	IO/LMCRootTreeWriter.cc
    	IO/LMCRootHistoWriter.cc
    	IO/LMCRootGraphWriter.cc
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
            ${PROJECT_SOURCE_DIR}/Source/Core/LMCRunParameters.hh 
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

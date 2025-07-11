    include_directories( BEFORE 
        ${CMAKE_CURRENT_SOURCE_DIR}
    )    

    set( testing_SOURCES
        testCatch.cc
        testFFT.cc      
        testMockEGun.cc
        testMockFreeField.cc
        testLMCTestSignal.cc
        testLMCPlaneWaveFIR.cc
        testLMCCavity.cc    
        testAliasingHF.cc    
    )
     
    set( lib_dependencies 
        LocustMC
    )     

    if( locust_mc_BUILD_WITH_ROOT ) 
        set( lib_dependencies ${lib_dependencies}
            LocustIO
        )
    endif( locust_mc_BUILD_WITH_ROOT )
    
     
    add_library( Catch INTERFACE )
    
    target_include_directories( Catch INTERFACE ${Scarab_LOCATION}/testing/catch )

    add_executable( RunTests ${testing_SOURCES} )
    target_link_libraries( RunTests ${lib_dependencies} )
    
    add_executable( testFFT testFFT.cc testCatch.cc)
    target_link_libraries( testFFT ${lib_dependencies})   
    
    add_executable( testMockEGun testMockEGun.cc testCatch.cc)
    target_link_libraries( testMockEGun ${lib_dependencies})    
    
    add_executable( testMockFreeField testMockFreeField.cc testCatch.cc)
    target_link_libraries( testMockFreeField ${lib_dependencies})

    add_executable( testLMCTestSignal testLMCTestSignal.cc testCatch.cc)
    target_link_libraries( testLMCTestSignal ${lib_dependencies})
        
    add_executable( testLMCPlaneWaveFIR testLMCPlaneWaveFIR.cc testCatch.cc)
    target_link_libraries( testLMCPlaneWaveFIR ${lib_dependencies})
    
    add_executable( testLMCCavity testLMCCavity.cc testCatch.cc)
    target_link_libraries( testLMCCavity ${lib_dependencies})  
            
    add_executable( testAliasingHF testAliasingHF.cc testCatch.cc)
    target_link_libraries( testAliasingHF ${lib_dependencies})    
            
    set( programs
        RunTests
    )    
       
    pbuilder_executables( programs lib_dependencies )
    install( TARGETS RunTests testFFT testMockEGun testMockFreeField testLMCTestSignal testLMCPlaneWaveFIR testLMCCavity testAliasingHF DESTINATION ${PROJECT_BINARY_DIRECTORY} )
    
    
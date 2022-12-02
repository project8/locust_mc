    set( CPP_DIR Cpp )
    
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
    )
     
    set( lib_dependencies 
        LocustMC
    )     
     
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
    
    
    set( CPP_DIR Cpp )
    
    include_directories( BEFORE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    )    

    set( testing_SOURCES
        testCatch.cc
        testFFT.cc      
        testMockEGun.cc
        testMockFreeField.cc
    )
     
    set( lib_dependencies 
        LocustMC
    )     
     
    add_library( Catch INTERFACE )
    target_include_directories( Catch INTERFACE ${Scarab_LOCATION}/testing/catch )

    add_executable( RunTests ${testing_SOURCES} )
    target_link_libraries( RunTests ${lib_dependencies} )


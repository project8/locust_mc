if( locust_mc_ENABLE_TESTING )

    set( programs
        testCatch
        testFFT
        testMockEGun
    )
    
    set( lib_dependencies 
        LocustMC
    )
    
    add_library( Catch INTERFACE )
    target_include_directories( Catch INTERFACE ${Scarab_LOCATION}/testing/catch )
        
    pbuilder_executables( programs lib_dependencies )

endif( locust_mc_ENABLE_TESTING )

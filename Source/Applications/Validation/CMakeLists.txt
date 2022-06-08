if( locust_mc_ENABLE_TESTING )

    set( lib_dependencies 
        LocustMC
    )
    
    set( programs
        test_new_app_with_options
        testFFT
        testMockEGun
        testMockFreeField
    )    
    
    pbuilder_executables( programs lib_dependencies )

endif( locust_mc_ENABLE_TESTING )

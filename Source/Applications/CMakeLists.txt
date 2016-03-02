if( locust_mc_ENABLE_EXECUTABLES )

    set( lib_dependencies 
        LocustMC
    )
    
    set( programs
 #       LocustSim
    )
    
    pbuilder_executables( programs lib_dependencies )

endif( locust_mc_ENABLE_EXECUTABLES )

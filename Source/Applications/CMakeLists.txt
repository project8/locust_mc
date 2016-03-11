if( locust_mc_ENABLE_EXECUTABLES )

    set( lib_dependencies 
        LocustMC
    )
    
    set( programs
        LocustSim
    )
    
    if( locust_mc_BUILD_WITH_KASSIOPEIA )
        set( programs ${programs}
            LMCKassiopeia
        )
    endif( locust_mc_BUILD_WITH_KASSIOPEIA )
    
    pbuilder_executables( programs lib_dependencies )

endif( locust_mc_ENABLE_EXECUTABLES )

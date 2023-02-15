    set( lib_dependencies 
        LocustMC        
    )
    
    if( locust_mc_BUILD_WITH_ROOT ) 
        set( lib_dependencies ${lib_dependencies}
            LocustIO
        )
    endif( locust_mc_BUILD_WITH_ROOT )
            
    
    set( programs
        LocustSim
    )
    
    if( locust_mc_BUILD_WITH_KASSIOPEIA )
        set( programs ${programs}
            LMCKassiopeia
        )
    endif( locust_mc_BUILD_WITH_KASSIOPEIA )
    
    pbuilder_executables( programs lib_dependencies )    
    
    if (locust_mc_ENABLE_TESTING)
        add_subdirectory( Testing )    
    endif (locust_mc_ENABLE_TESTING)
    

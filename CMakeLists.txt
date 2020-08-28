cmake_minimum_required( VERSION 3.1 )

#########
# setup #
#########

# Define the project
cmake_policy( SET CMP0048 NEW ) # version in project()
project( locust_mc VERSION 2.1.2)

list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/Scarab/cmake )
include( PackageBuilder )

pbuilder_prepare_project()

set_option( USE_CPP11 TRUE )

option( locust_mc_BUILD_WITH_KASSIOPEIA "Option to build with Kassiopeia" FALSE )
option( locust_mc_BUILD_WITH_ROOT "Option to build with ROOT" TRUE )

# make the data install directory available as a preprocessor macro
add_definitions( -DPB_DATA_INSTALL_DIR=${DATA_INSTALL_DIR} )

set (LOCUST_MC_INCLUDE_DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Transmitters/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/RxComponents/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Transforms/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/IO/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Distributions/
)


if (locust_mc_BUILD_WITH_KASSIOPEIA)

    # Kassiopeia needs to be built independently from, but as a submodule and dependency of, Locust.  
    # In other words, the configuration and build environments need to be kept separate.
    # However, Locust needs to build against Kassiopeia.
    # When Kassiopeia is run from within Locust it knows about the Kassiopeia-derived Locust classes.

    # First we run the configure stage of Kassiopeia during the configure (cmake) process of Locust.
    # This build takes place in its own build directory, but it installs in the same place as Locust.
    # We run cmake twice so that we get all of the variables setup correctly.
    # We do this here instead of below in the ExternalProject so that Locust can learn about the Kassiopeia targets at the configure stage.
    execute_process( 
        COMMAND mkdir kassiopeia
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
    execute_process( 
        COMMAND cmake -D CMAKE_INSTALL_PREFIX:STRING=${CMAKE_INSTALL_PREFIX} 
                      -D CMAKE_INSTALL_LIBDIR:STRING=lib
                      -D BUILD_KASSIOPEIA:BOOL=TRUE 
                      -D BUILD_KEMFIELD:BOOL=TRUE 
                      -D BUILD_KGEOBAG:BOOL=TRUE 
                      -D BUILD_KOMMON:BOOL=TRUE 
                      ${PROJECT_SOURCE_DIR}/kassiopeia
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/kassiopeia
    )
    #execute_process( 
    #    COMMAND cmake -D CMAKE_INSTALL_PREFIX:STRING=${CMAKE_INSTALL_PREFIX} -D BUILD_KASSIOPEIA:BOOL=TRUE -D BUILD_KEMFIELD:BOOL=TRUE -D BUILD_KGEOBAG:BOOL=TRUE -D BUILD_KOMMON:BOOL=TRUE ${PROJECT_SOURCE_DIR}/kassiopeia
    #    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/kassiopeia
    #)

    # Second we use the ExternalProject module to build and install Kassiopeia.
    # This takes place during the build phase of Locust.
    # Note that we disable the configure stage of the ExternalProject because we already took care of the configure stage above.
    include(ExternalProject)
    
    ExternalProject_Add( KassiopeiaExt
        PREFIX            ${CMAKE_INSTALL_PREFIX}
        SOURCE_DIR        ${PROJECT_SOURCE_DIR}/kassiopeia
        BINARY_DIR        ${PROJECT_BINARY_DIR}/kassiopeia
        CONFIGURE_COMMAND ""
        BUILD_ALWAYS      1
        BUILD_COMMAND     make
        INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
        INSTALL_COMMAND   make -j3 install
    )

    # Finally we use find_package() to learn about the Kassiopeia build targets.
    find_package( Kommon REQUIRED CONFIG HINTS ${PROJECT_BINARY_DIR}/kassiopeia/Kommon )
    find_package( KGeoBag REQUIRED CONFIG HINTS ${PROJECT_BINARY_DIR}/kassiopeia/KGeoBag )
    find_package( KEMField REQUIRED CONFIG HINTS ${PROJECT_BINARY_DIR}/kassiopeia/KEMField )
    find_package( Kassiopeia REQUIRED CONFIG HINTS ${PROJECT_BINARY_DIR}/kassiopeia/Kassiopeia )

endif (locust_mc_BUILD_WITH_KASSIOPEIA)


#####################
# locust submodules #
#####################

pbuilder_add_submodule( Scarab Scarab/library )

pbuilder_add_submodule( Monarch monarch )
# we need the definitions for HAS_ATTR_IFC and HAS_GRP_IFC from Monarch to use the proper class interfaces in M3Header and M3Stream
set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} APPEND PROPERTY COMPILE_DEFINITIONS ${Monarch_COMPILE_DEFINITIONS} )


#######################
# locust dependencies #
#######################

find_package(GSL REQUIRED)
if (GSL_FOUND)
    pbuilder_add_ext_libraries( ${GSL_LIBRARIES} )
else (GSL_FOUND)
    message(FATAL "GSL not found")
endif (GSL_FOUND)


# FFTW
find_package(FFTW REQUIRED)
if (NOT FFTW_THREADS_FOUND)
  message(FATAL "FFTW_THREADS not found")
endif (NOT FFTW_THREADS_FOUND)
    pbuilder_add_ext_libraries (${FFTW_LIBRARIES})
    set (FFTW_NTHREADS 1 CACHE STRING "Number of threads to use for FFTW processes")
    add_definitions (-DFFTW_NTHREADS=${FFTW_NTHREADS})
    message (STATUS "FFTW configured to use up to ${FFTW_NTHREADS} threads.")
message (STATUS "FFTW_LIBRARIES is ${FFTW_LIBRARIES}")
include_directories (${FFTW_INCLUDE_DIR})

# ROOT
if (locust_mc_BUILD_WITH_ROOT)
    find_package (ROOT 6.00 COMPONENTS Gui Spectrum TMVA)
else (locust_mc_BUILD_WITH_ROOT)
    set (ROOT_FOUND FALSE)
endif (locust_mc_BUILD_WITH_ROOT)
if (ROOT_FOUND)
    add_definitions(-DROOT_FOUND)
    pbuilder_add_ext_libraries (${ROOT_LIBRARIES})
    message(STATUS "ROOT_LIBRARIES is ${ROOT_LIBRARIES}")
else (ROOT_FOUND)
    message(STATUS "Building without ROOT")
    remove_definitions(-DROOT_FOUND)
endif (ROOT_FOUND)
include_directories (${ROOT_INCLUDE_DIR})


###########
# library #
###########

if( locust_mc_BUILD_WITH_KASSIOPEIA )
    # Here we add the include directories, link directories, and libraries that will be used by Kassiopeia to Locust's build configuration.

    set (LOCUST_MC_INCLUDE_DIRECTORIES
        ${LOCUST_MC_INCLUDE_DIRECTORIES}
        ${CMAKE_CURRENT_SOURCE_DIR}/Source/Kassiopeia/
        ${Kommon_INCLUDE_DIRS}
        ${KGeoBag_INCLUDE_DIRS}
        ${KEMField_INCLUDE_DIRS}
        ${Kassiopeia_INCLUDE_DIRS}
    )

    set (LOCUST_MC_LINK_DIRECTORIES
        ${LOCUST_MC_LINK_DIRECTORIES}
        ${LIB_INSTALL_DIR}        
    )

    # We need to have the LocustMC library depend on the Kassiopeia libraries, and to the KassiopeiaExt target
    # To do this we create an imported library for each Kassiopeia library
    # We then add a dependency of KassiopeiaExt to each imported library target
    # LocustMC library will link to those imported libraries, and therefore depend on KassiopeiaExt

    set( All_Kassiopeia_LIBRARIES ${Kommon_LIBRARIES} ${KGeoBag_LIBRARIES} ${KEMField_LIBRARIES} ${Kassiopeia_LIBRARIES} )
    foreach( LIBRARY IN LISTS All_Kassiopeia_LIBRARIES )
        add_library( ${LIBRARY} SHARED IMPORTED )
        #message( STATUS "****** constructed name: ${LIB_INSTALL_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}${LIBRARY}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
        set_target_properties( ${LIBRARY} PROPERTIES IMPORTED_LOCATION "${LIB_INSTALL_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}${LIBRARY}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
        add_dependencies( ${LIBRARY} KassiopeiaExt )
    endforeach()
      
              
    pbuilder_add_ext_libraries(
        ${Kommon_LIBRARIES}
        ${KGeoBag_LIBRARIES}
        ${KEMField_LIBRARIES}
        ${Kassiopeia_LIBRARIES}
    )

endif( locust_mc_BUILD_WITH_KASSIOPEIA )

include_directories( ${LOCUST_MC_INCLUDE_DIRECTORIES} )
link_directories( ${LOCUST_MC_LINK_DIRECTORIES} )

add_subdirectory( Source )
add_subdirectory( Data )
add_subdirectory( Data/TransferFunctions )

###############
# executables #
###############

if (locust_mc_ENABLE_EXECUTABLES )
    add_subdirectory( Source/Applications )
endif (locust_mc_ENABLE_EXECUTABLES )

################
# Config files #
################

add_subdirectory( Config )


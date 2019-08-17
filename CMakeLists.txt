cmake_minimum_required( VERSION 3.1 )

#########
# setup #
#########

# Define the project
cmake_policy( SET CMP0048 NEW ) # version in project()
project( locust_mc VERSION 1.13.3)

list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/Scarab/cmake )
include( PackageBuilder )

pbuilder_prepare_project()


set_option( USE_CPP11 TRUE )

option( locust_mc_BUILD_WITH_KASSIOPEIA "Option to build with Kassiopeia" FALSE )
option( locust_mc_BUILD_WITH_ROOT "Option to build with ROOT" TRUE )


set (LOCUST_MC_INCLUDE_DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Transmitters/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Transforms/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/
)


if (locust_mc_BUILD_WITH_KASSIOPEIA)

    include(ExternalProject)
    message( STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}" )
    
    ExternalProject_Add( KassiopeiaExt
        PREFIX            ${CMAKE_INSTALL_PREFIX}
        SOURCE_DIR        ${PROJECT_SOURCE_DIR}/kassiopeia
        CMAKE_ARGS        -U *_INSTALL_DIR
        CMAKE_CACHE_ARGS  -D CMAKE_INSTALL_PREFIX:STRING=${CMAKE_INSTALL_PREFIX} -D BUILD_KASSIOPEIA:BOOL=TRUE -D BUILD_KEMFIELD:BOOL=TRUE -D BUILD_KGEOBAG:BOOL=TRUE -D BUILD_KOMMON:BOOL=TRUE
        BINARY_DIR        ${PROJECT_BINARY_DIR}/kassiopeia
        CMAKE_COMMAND     cmake
        BUILD_ALWAYS      1
        BUILD_COMMAND     make
        INSTALL_DIR       ${CMAKE_INSTALL_PREFIX}
        INSTALL_COMMAND   make install
    )

    ExternalProject_Add_Step( KassiopeiaExt prelim_configure
        COMMAND           cmake -D CMAKE_INSTALL_PREFIX:STRING=${CMAKE_INSTALL_PREFIX} -D BUILD_KASSIOPEIA:BOOL=TRUE -D BUILD_KEMFIELD:BOOL=TRUE -D BUILD_KGEOBAG:BOOL=TRUE -D BUILD_KOMMON:BOOL=TRUE ${PROJECT_SOURCE_DIR}/kassiopeia
        DEPENDERS         configure
        BYPRODUCTS        ${PROJECT_BINARY_DIR}/kassiopeia/CMakeCache.txt
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/kassiopeia
        USES_TERMINAL     1
    )

    set (LOCUST_MC_INCLUDE_DIRECTORIES
        ${LOCUST_MC_INCLUDE_DIRECTORIES}
        ${CMAKE_CURRENT_SOURCE_DIR}/Source/Kassiopeia/
        ${CMAKE_INSTALL_PREFIX}/include
        ${INCLUDE_INSTALL_DIR}
    )

    set (LOCUST_MC_LINK_DIRECTORIES
        ${LIB_INSTALL_DIR}
    )

    find_package( Kommon REQUIRED CONFIG HINTS ${LIB_INSTALL_DIR}/cmake/Kommon )
    find_package( KGeoBag REQUIRED CONFIG HINTS ${LIB_INSTALL_DIR}/cmake/KGeoBag )
    find_package( KEMField REQUIRED CONFIG HINTS ${LIB_INSTALL_DIR}/cmake/KEMField )
    find_package( Kassiopeia REQUIRED CONFIG HINTS ${LIB_INSTALL_DIR}/cmake/Kassiopeia )

    pbuilder_add_ext_libraries(
        ${Kommon_LIBRARIES}
        ${KGeoBag_LIBRARIES}
        ${KEMField_LIBRARIES}
        ${Kassiopeia_LIBRARIES}
    )

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

###############
# executables #
###############
if (locust_mc_ENABLE_EXECUTABLES )
    add_subdirectory( Source/Applications )
endif (locust_mc_ENABLE_EXECUTABLES )






#pbuilder_install_headers(${CMAKE_CURRENT_SOURCE_DIR}/Source/ReceiverTransferFunctions.hh)

#add_executable( generate_data ${CMAKE_CURRENT_SOURCE_DIR}/Source/generate_data.cc  ${CMAKE_CURRENT_SOURCE_DIR}/Source/ReceiverTransferFunctions.cc )
#target_link_libraries( generate_data ${EXTERNAL_LIBRARIES} )
#pbuilder_install_executables( generate_data )

#add_executable( check_egg_digitization ${CMAKE_CURRENT_SOURCE_DIR}/Source/check_egg_digitization.cc )
#target_link_libraries( check_egg_digitization ${EXTERNAL_LIBRARIES} )
#pbuilder_install_executables( check_egg_digitization )

pbuilder_install_config (

)


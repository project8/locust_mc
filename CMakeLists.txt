cmake_minimum_required( VERSION 3.1 )

#########
# setup #
#########

# Define the project
cmake_policy( SET CMP0048 NEW ) # version in project()
project( locust_mc VERSION 1.15.9)

list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/Scarab/cmake )
include( PackageBuilder )

pbuilder_prepare_project()


set_option( USE_CPP11 TRUE )

option( locust_mc_BUILD_WITH_KASSIOPEIA "Option to build with Kassiopeia" FALSE )
option( locust_mc_BUILD_WITH_ROOT "Option to build with ROOT" TRUE )

# make the data install directory available as a preprocessor macro
add_definitions( -DPB_DATA_INSTALL_DIR=${DATA_INSTALL_DIR} )


#####################
# locust submodules #
#####################


pbuilder_add_submodule( Scarab Scarab/library )

pbuilder_add_submodule( Monarch monarch )
# we need the definitions for HAS_ATTR_IFC and HAS_GRP_IFC from Monarch to use the proper class interfaces in M3Header and M3Stream
set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} APPEND PROPERTY COMPILE_DEFINITIONS ${Monarch_COMPILE_DEFINITIONS} )

# set this up before the Kassiopeia build, so that the Kassiopeia build can use Locust include directories
set (LOCUST_MC_INCLUDE_DIRECTORIES
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Transmitters/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/RxComponents/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Transforms/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Kassiopeia/
)

include_directories( BEFORE ${LOCUST_MC_INCLUDE_DIRECTORIES} )


if (locust_mc_BUILD_WITH_KASSIOPEIA)

    # Keep KasperDefaults from overriding the install prefix
    set( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE )

    list( APPEND CMAKE_MODULE_PATH 
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kommon/cmake 
    )
    include( KasperDefaults )
    set( Kommon_DIR ${PROJECT_BINARY_DIR}/kassiopeia/Kommon )
    set( KGeoBag_DIR ${PROJECT_BINARY_DIR}/kassiopeia/KGeoBag )
    set( KEMField_DIR ${PROJECT_BINARY_DIR}/kassiopeia/KEMField )

    set( CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)    
    set( DATA_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/data )    

    add_subdirectory( kassiopeia )
    
    if( Kommon_ENABLE_DEBUG )
        add_definitions( -DKommon_ENABLE_DEBUG )
    endif( Kommon_ENABLE_DEBUG )
    if( Kommon_USE_VTK )
        add_definitions( -DKommon_ENABLE_DEBUG )
    endif( Kommon_USE_VTK )
    if( Kommon_USE_ROOT )
        add_definitions( -DKommon_USE_ROOT )
    endif( Kommon_USE_ROOT )
    
    if( KGeoBag_ENABLE_DEBUG )
        add_definitions( -DKGeoBag_ENABLE_DEBUG )
    endif( KGeoBag_ENABLE_DEBUG )
    if( KGeoBag_USE_ROOT )
        add_definitions( -DKGeoBag_USE_ROOT )
    endif( KGeoBag_USE_ROOT )
    if( KGeoBag_USE_VTK )
        add_definitions( -DKGeoBag_USE_VTK )
    endif( KGeoBag_USE_VTK )
    
    if( KEMField_USE_ROOT )
        add_definitions( -DKEMField_USE_ROOT )
    endif( KEMField_USE_ROOT )
    if( KEMField_USE_VTK )
        add_definitions( -DKEMField_USE_VTK )
    endif( KEMField_USE_VTK )
    
    if( Kassiopeia_ENABLE_DEBUG )
        add_definitions( -DKassiopeia_ENABLE_DEBUG )
    endif( Kassiopeia_ENABLE_DEBUG )
    if( Kassiopeia_USE_ROOT )
        add_definitions( -DKassiopeia_USE_ROOT )
    endif( Kassiopeia_USE_ROOT )
    if( Kassiopeia_USE_VTK )
        add_definitions( -DKassiopeia_USE_VTK )
    endif( Kassiopeia_USE_VTK )
    
    include( ${PROJECT_BINARY_DIR}/kassiopeia/Kommon/KommonConfig.cmake )
    include_directories( ${Kommon_INCLUDE_DIRS} )
    pbuilder_add_ext_libraries( ${Kommon_LIBRARIES} )
    
    include( ${PROJECT_BINARY_DIR}/kassiopeia/KGeoBag/KGeoBagConfig.cmake )
    include_directories( ${KGeoBag_INCLUDE_DIRS} )
    pbuilder_add_ext_libraries( ${KGeoBag_LIBRARIES} )
    
    include( ${PROJECT_BINARY_DIR}/kassiopeia/KEMField/KEMFieldConfig.cmake )
    include_directories( ${KEMField_INCLUDE_DIRS} )
    pbuilder_add_ext_libraries( ${KEMField_LIBRARIES} )
    
    include( ${PROJECT_BINARY_DIR}/kassiopeia/Kassiopeia/KassiopeiaConfig.cmake )
    include_directories( ${Kassiopeia_INCLUDE_DIRS} )
    pbuilder_add_ext_libraries( ${Kassiopeia_LIBRARIES} )
endif (locust_mc_BUILD_WITH_KASSIOPEIA)


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

add_subdirectory( Source )
add_subdirectory( Data )



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


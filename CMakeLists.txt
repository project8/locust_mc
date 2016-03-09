cmake_minimum_required( VERSION 3.0 )

#########
# setup #
#########

# Define the project
cmake_policy( SET CMP0048 NEW ) # version in project()
project( locust_mc VERSION 1.0.0 )

list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/Scarab/cmake )
include( PackageBuilder )

pbuilder_prepare_project()


set_option( USE_CPP11 FALSE )
set( CMAKE_CXX_FLAGS "-std=c++11")  # workaround in case the above line does not work.

option( locust_mc_BUILD_WITH_KASSIOPEIA "Option to build with Kassiopeia" FALSE )

# cxx flags to avoid incompatibility between gcc 5 and gcc <5.
#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_INIT} -Wall -Wextra -D_GLIBCXX_USE_CXX11_ABI=0")


list ( APPEND CMAKE_MODULE_PATH 
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/kassiopeia/Kommon/cmake 
)


if (locust_mc_BUILD_WITH_KASSIOPEIA)
    include( KasperDefaults )
    set( KASPER_BUILD TRUE )
endif (locust_mc_BUILD_WITH_KASSIOPEIA)

set ( CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR} )


list ( APPEND CMAKE_PREFIX_PATH 
    ${CMAKE_BINARY_DIR}/kassiopeia/Kommon
    ${CMAKE_BINARY_DIR}/kassiopeia/KGeoBag
    ${CMAKE_BINARY_DIR}/kassiopeia/KEMField
)



IF (locust_mc_BUILD_WITH_KASSIOPEIA)

    # use, i.e. don't skip the full RPATH for the build tree
    set( CMAKE_SKIP_BUILD_RPATH TRUE )
    
    # when building, don't use the install RPATH already
    #( but later on when installing )
    set( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
    
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
    
    set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
    
    # the RPATH to be used when installing, but only if it's not a system directory
    LIST( FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir )
    IF( "${isSystemDir}" STREQUAL "-1" )
       set( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
    ENDIF( "${isSystemDir}" STREQUAL "-1" )
    
    
    
    
    
    option( GCC_FIX_LINKER_ERROR "Fix linker errors with some GCC versions by adding the --no-as-needed flag" ON )
    if( GCC_FIX_LINKER_ERROR )
        set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed" )
    endif( GCC_FIX_LINKER_ERROR )



    set (KASSIOPEIA_INCLUDE_DIRECTORIES
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/ 
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Operators/Include/ 
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/ 
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kommon/Core/Utility/ 
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kommon/Core/Initialization/ 
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kommon/Core/File/ 
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kommon/Extensions/Root/File/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kommon/Extensions/Root/Initialization/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KGeoBag/Source/Math/Include/ 
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Utility/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kommon/Extensions/Gsl/Utility/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Fields/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Generators/Include/ 
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Geometry/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Interactions/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Math/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Navigators/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Modifiers/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Terminators/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Trajectories/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Visualization/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Writers/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Readers/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/Plugins/KGeoBag/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/Surfaces/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/Core/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/Math/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/Math/Array/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/Math/Utilities/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/IO/Streamers/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/IO/HashGenerator/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/IO/FileManipulation/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/IO/StructuredASCII/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/BoundaryIntegrals/Electrostatic/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/LinearAlgebra/Core/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/LinearAlgebra/Solvers/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/LinearAlgebra/Visitors/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/BoundaryIntegrals/Core/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FieldSolvers/Integrating/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FieldSolvers/ZonalHarmonic/Solver/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FieldSolvers/ZonalHarmonic/Generator/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/ExternalFields/Electromagnets/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FastMultipole/Interface/BoundaryIntegrals/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FastMultipole/Electrostatics/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FastMultipole/Core/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FastMultipole/Tree/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FastMultipole/Math/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FastMultipole/Kernel/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FastMultipole/Interface/FieldSolvers/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KEMField/Source/2.0/FastMultipole/Interface/Extraction/include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KGeoBag/Source/Core/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KGeoBag/Source/Shapes/Complex/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KGeoBag/Source/Shapes/Utility/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KGeoBag/Source/Shapes/RotatedVolumes/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KGeoBag/Source/Shapes/PlanarShapes/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KGeoBag/Source/Shapes/RotatedAreas/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KGeoBag/Source/Extensions/Mesh/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KGeoBag/Source/Extensions/AxialMesh/Include/
        ${PROJECT_SOURCE_DIR}/kassiopeia/KGeoBag/Source/Extensions/DiscreteRotationalMesh/Include/
        /usr/local/root/include/
    )
    include_directories( ${KASSIOPEIA_INCLUDE_DIRECTORIES} )

ENDIF (locust_mc_BUILD_WITH_KASSIOPEIA)


set (LOCUST_MC_INCLUDE_DIRECTORIES
    ${PROJECT_SOURCE_DIR}/Source/Generators/
    ${PROJECT_SOURCE_DIR}/Source/Core/
)
include_directories( ${LOCUST_MC_INCLUDE_DIRECTORIES} )



pbuilder_add_submodule( Scarab Scarab/library )

pbuilder_add_submodule( Monarch monarch )





#######################
# locust dependencies #
#######################



include_directories( ${PROJECT_SOURCE_DIR}/Scarab/library/param/codec/json/RapidJSON )

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



IF (locust_mc_BUILD_WITH_KASSIOPEIA)
set( CMAKE_CXX_FLAGS " ")
    add_subdirectory( kassiopeia/Kommon )
    #include_directories( ${Kommon_INCLUDE_DIRS} )
    #link_directories( ${Kommon_LIBRARY_DIRS} )
    file(GLOB_RECURSE Kommon_LIBRARIES ${CMAKE_BINARY_DIR}/kassiopeia/Kommon/*.so)
#    pbuilder_add_ext_libraries( ${Kommon_LIBRARIES} )
#    pbuilder_add_ext_libraries( ${KommonCore} )
    message( STATUS "Kommon libraries: ${Kommon_LIBRARIES}" )    
    
    set( BUILD_KOMMON ON CACHE BOOL "${BUILD_KOMMON_DESC} (Required)" FORCE )
    set( Kommon_USE_GSL ON CACHE BOOL "(Required)" FORCE )
    set( Kommon_USE_ROOT ON CACHE BOOL "(Required)" FORCE )
    set( Kommon_USE_VTK ON CACHE BOOL "(Required)" FORCE )
    add_subdirectory( kassiopeia/KGeoBag )
    include_directories( ${KGeoBag_INCLUDE_DIRS} )
    link_directories( ${KGeoBag_LIBRARY_DIRS} )
    file(GLOB_RECURSE KGeoBag_LIBRARIES ${CMAKE_BINARY_DIR}/kassiopeia/KGeoBag/*.so)
#    pbuilder_add_ext_libraries( ${KGeoBag_LIBRARIES} )
#    pbuilder_add_ext_libraries( ${KGeoBagCore} )
    message( STATUS "KGeoBag libraries: ${KGeoBag_LIBRARIES}" )
    

    add_subdirectory( kassiopeia/KEMField )
    include_directories( ${KEMField_INCLUDE_DIRS} )
    link_directories( ${KEMField_LIBRARY_DIRS} )
    file(GLOB_RECURSE KEMField_LIBRARIES ${CMAKE_BINARY_DIR}/kassiopeia/KEMField/*.so)
#    pbuilder_add_ext_libraries( ${KEMField_LIBRARIES} )
#    pbuilder_add_ext_libraries( ${KEMCore} )
    message( STATUS "KEMField libraries: ${KEMField_LIBRARIES}" )

    
    set( BUILD_KOMMON ON CACHE BOOL "${BUILD_KOMMON_DESC} (Required)" FORCE )
    set( BUILD_KGEOBAG ON CACHE BOOL "${BUILD_KGEOBAG_DESC} (Required)" FORCE )
    set( BUILD_KEMFIELD ON CACHE BOOL "${BUILD_KEMFIELD_DESC} (Required)" FORCE )
    set( Kommon_USE_GSL ON CACHE BOOL "(Required)" FORCE )
    set( Kommon_USE_ROOT ON CACHE BOOL "(Required)" FORCE )
    set( KGeoBag_USE_ROOT ON CACHE BOOL "(Required)" FORCE )
    set( Kommon_USE_VTK ON CACHE BOOL "(Required)" FORCE )
    set( KGeoBag_USE_VTK ON CACHE BOOL "(Required)" FORCE )
    add_subdirectory( kassiopeia/Kassiopeia )
    include_directories( ${Kassiopeia_INCLUDE_DIRS} )
    link_directories( ${Kassiopeia_LIBRARY_DIRS} )
    file(GLOB_RECURSE Kassiopeia_LIBRARIES ${CMAKE_BINARY_DIR}/kassiopeia/Kassiopeia/*.so)
#    pbuilder_add_ext_libraries( ${Kassiopeia_LIBRARIES} )
    message( STATUS "Kassiopeia_LIBRARIES: ${Kassiopeia_LIBRARIES}" )

set( CMAKE_CXX_FLAGS "-std=c++11") # switch cpp11 support back on.
ENDIF ( locust_mc_BUILD_WITH_KASSIOPEIA )

message( STATUS "cmake_binary_dir is  ${CMAKE_BINARY_DIR}" )
    message(STATUS "cmake_exe_linker_flags:  ${CMAKE_EXE_LINKER_FLAGS}" )

#pbuilder_add_ext_libraries( "m" )



###########
# library #
###########

add_subdirectory( Source )

###############
# executables #
###############

add_subdirectory( Source/Applications )






#pbuilder_install_headers(${CMAKE_CURRENT_SOURCE_DIR}/Source/ReceiverTransferFunctions.hh)

#add_executable( generate_data ${CMAKE_CURRENT_SOURCE_DIR}/Source/generate_data.cc  ${CMAKE_CURRENT_SOURCE_DIR}/Source/ReceiverTransferFunctions.cc )
#target_link_libraries( generate_data ${EXTERNAL_LIBRARIES} )
#pbuilder_install_executables( generate_data )

#add_executable( check_egg_digitization ${CMAKE_CURRENT_SOURCE_DIR}/Source/check_egg_digitization.cc )
#target_link_libraries( check_egg_digitization ${EXTERNAL_LIBRARIES} )
#pbuilder_install_executables( check_egg_digitization )

pbuilder_install_config (
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test.json          
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/receiver_transfer_functions.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_bg_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_bg_single_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_single_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/tone_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/tone_single_config.json
)


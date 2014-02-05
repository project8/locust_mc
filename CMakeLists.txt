cmake_minimum_required( VERSION 2.8 )

#########
# setup #
#########

project( locust_mc )

list (APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

include( ${PROJECT_SOURCE_DIR}/PackageBuilder.cmake )
pbuilder_prepare_project( 1 0 0 )
include_directories( ${PROJECT_SOURCE_DIR} )


#######################
# locust dependencies #
#######################

pbuilder_add_ext_libraries( "m" )

include_directories( ${PROJECT_SOURCE_DIR}/RapidJSON )

add_subdirectory( RandomLib )
include_directories( ${PROJECT_SOURCE_DIR}/RandomLib/include )

#include_directories( ${PROJECT_SOURCE_DIR}/MTwist )
#set( MTWIST_SOURCE_FILES
#    ${PROJECT_SOURCE_DIR}/MTwist/mtwist.c
#    ${PROJECT_SOURCE_DIR}/MTwist/randistrs.c
#)

add_subdirectory( monarch )
pbuilder_add_ext_libraries( ${Monarch_LIBRARIES} )
# add to the RPATH to be used when installing, but only if it's not a system directory
list( FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${Monarch_LIBRARY_DIR}" isSystemDir )
if( "${isSystemDir}" STREQUAL "-1" )
   list( APPEND CMAKE_INSTALL_RPATH ${Monarch_LIBRARY_DIR} )
endif( "${isSystemDir}" STREQUAL "-1" )
include_directories( ${PROJECT_SOURCE_DIR}/monarch/Include )
include_directories( ${PROJECT_SOURCE_DIR}/monarch/libthorax/cpp )
include_directories( ${PROJECT_SOURCE_DIR}/monarch/libthorax )

#find_package( Yajl )
#if (NOT YAJL_FOUND)
#  message(FATAL "Yajl not found")
#endif (NOT YAJL_FOUND)
#include_directories( ${YAJL_INCLUDE_DIR} )
#pbuilder_add_ext_libraries( ${YAJL_LIBRARIES} )

find_package( Protobuf )
include_directories( ${PROTOBUF_INCLUDE_DIR} )
pbuilder_add_ext_libraries( ${PROTOBUF_LIBRARIES} )


# FFTW
find_package(FFTW REQUIRED)
if (NOT FFTW_THREADS_FOUND)
  message(FATAL "FFTW_THREADS not found")
endif (NOT FFTW_THREADS_FOUND)
    pbuilder_add_ext_libraries (${FFTW_LIBRARIES})
    set (FFTW_NTHREADS 1 CACHE STRING "Number of threads to use for FFTW processes")
    add_definitions (-DFFTW_NTHREADS=${FFTW_NTHREADS})
    message (STATUS "FFTW configured to use up to ${FFTW_NTHREADS} threads.")
include_directories (${FFTW_INCLUDE_DIR})


###########
# library #
###########

set( LOCUST_MC_HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCConfigurator.hh
#    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCConstants.hh
#    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCDefaultConfig.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCException.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCGaussianNoiseGenerator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCGenerator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCLogger.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCParam.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCParser.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCSignal.hh
    #${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCSingleAmpSim.hh
)

set( LOCUST_MC_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCConfigurator.cc
#    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCDefaultConfig.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCException.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCGaussianNoiseGenerator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCGenerator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCLogger.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCParam.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCParser.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCSignal.cc
    #${CMAKE_CURRENT_SOURCE_DIR}/Source/LMCSingleAmpSim.cc
)

add_library( LocustMC ${LOCUST_MC_SOURCE_FILES} )
target_link_libraries( LocustMC ${EXTERNAL_LIBRARIES} )

pbuilder_install_headers( ${LOCUST_MC_HEADER_FILES} )
pbuilder_install_libraries( LocustMC )


###############
# executables #
###############

#pbuilder_install_headers(${CMAKE_CURRENT_SOURCE_DIR}/Source/ReceiverTransferFunctions.hh)

#add_executable( generate_data ${CMAKE_CURRENT_SOURCE_DIR}/Source/generate_data.cc  ${CMAKE_CURRENT_SOURCE_DIR}/Source/ReceiverTransferFunctions.cc )
#target_link_libraries( generate_data ${EXTERNAL_LIBRARIES} )
#pbuilder_install_executables( generate_data )

#add_executable( check_egg_digitization ${CMAKE_CURRENT_SOURCE_DIR}/Source/check_egg_digitization.cc )
#target_link_libraries( check_egg_digitization ${EXTERNAL_LIBRARIES} )
#pbuilder_install_executables( check_egg_digitization )

pbuilder_install_config (
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/receiver_transfer_functions.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_bg_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_bg_single_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_single_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/tone_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/tone_single_config.json
)
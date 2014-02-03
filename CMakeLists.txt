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

find_package( Yajl )
if (NOT YAJL_FOUND)
  message(FATAL "Yajl not found")
endif (NOT YAJL_FOUND)
include_directories( ${YAJL_INCLUDE_DIR} )
pbuilder_add_ext_libraries( ${YAJL_LIBRARIES} )

find_package( Protobuf )
include_directories( ${PROTOBUF_INCLUDE_DIR} )
pbuilder_add_ext_libraries( ${PROTOBUF_LIBRARIES} )


# FFTW
find_package(FFTW REQUIRED)
if (NOT FFTW_THREADS_FOUND)
  message(FATAL "FFTW_THREADS not found")
endif (NOT FFTW_THREADS_FOUND)
    pbuilder_add_ext_libraries (${FFTW_LIBRARIES})
    if (FFTW_THREADS_FOUND AND NOT Katydid_SINGLETHREADED)
        set (FFTW_NTHREADS 1 CACHE STRING "Number of threads to use for FFTW processes")
        add_definitions (-DFFTW_NTHREADS=${FFTW_NTHREADS})
        message (STATUS "FFTW configured to use up to ${FFTW_NTHREADS} threads.")
    else (FFTW_THREADS_FOUND AND NOT Katydid_SINGLETHREADED)
        remove_definitions (-DFFTW_NTHREADS=${FFTW_NTHREADS})
    endif (FFTW_THREADS_FOUND AND NOT Katydid_SINGLETHREADED)
include_directories (${FFTW_INCLUDE_DIR})


###############
# executables #
###############

pbuilder_install_headers(${CMAKE_CURRENT_SOURCE_DIR}/Source/ReceiverTransferFunctions.hh)

add_executable( generate_data ${CMAKE_CURRENT_SOURCE_DIR}/Source/generate_data.cc  ${CMAKE_CURRENT_SOURCE_DIR}/Source/ReceiverTransferFunctions.cc )
target_link_libraries( generate_data ${EXTERNAL_LIBRARIES} )
pbuilder_install_executables( generate_data )

add_executable( check_egg_digitization ${CMAKE_CURRENT_SOURCE_DIR}/Source/check_egg_digitization.cc )
target_link_libraries( check_egg_digitization ${EXTERNAL_LIBRARIES} )
pbuilder_install_executables( check_egg_digitization )

pbuilder_install_config (
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/receiver_transfer_functions.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_bg_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_bg_single_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/test_single_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/tone_config.json
    ${CMAKE_CURRENT_SOURCE_DIR}/Config/tone_single_config.json
)
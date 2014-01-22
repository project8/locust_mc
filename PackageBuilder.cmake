# PackageBuilder.cmake
# Author: Noah Oblath
# Parts of this script are based on work done by Sebastian Voecking and Marco Haag in the Kasper package
# Convenient macros and default variable settings for the Katydid build.

# check if this is a stand-alone build
set( PBUILDER_STANDALONE FALSE CACHE INTERNAL "Flag for whether or not this is a stand-alone build" )
if( ${CMAKE_SOURCE_DIR} STREQUAL ${PROJECT_SOURCE_DIR} )
    set( PBUILDER_STANDALONE TRUE )
endif( ${CMAKE_SOURCE_DIR} STREQUAL ${PROJECT_SOURCE_DIR} )

# preprocessor defintion for debug build
if ("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
    add_definitions(-D${PROJECT_NAME}_DEBUG)
else ("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
    remove_definitions(-D${PROJECT_NAME}_DEBUG)    
endif ("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")

# Setup the default install prefix
# This gets set to the binary directory upon first configuring.
# If the user changes the prefix, but leaves the flag OFF, then it will remain as the user specified.
# If the user wants to reset the prefix to the default (i.e. the binary directory), then the flag should be set ON.
if (NOT DEFINED SET_INSTALL_PREFIX_TO_DEFAULT)
    set (SET_INSTALL_PREFIX_TO_DEFAULT ON)
endif (NOT DEFINED SET_INSTALL_PREFIX_TO_DEFAULT)
if (SET_INSTALL_PREFIX_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE PATH "Install prefix" FORCE)
    set (SET_INSTALL_PREFIX_TO_DEFAULT OFF CACHE BOOL "Reset default install path when when configuring" FORCE)
endif (SET_INSTALL_PREFIX_TO_DEFAULT)

# install subdirectories
set (INCLUDE_INSTALL_SUBDIR "include" CACHE PATH "Install subdirectory for headers")
set (LIB_INSTALL_SUBDIR "lib" CACHE PATH "Install subdirectory for libraries")
set (BIN_INSTALL_SUBDIR "bin" CACHE PATH "Install subdirectory for binaries")
set (CONFIG_INSTALL_SUBDIR "config" CACHE PATH "Install subdirectory for config files")

set (INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_SUBDIR}")
set (LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_SUBDIR}")
set (BIN_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${BIN_INSTALL_SUBDIR}")
set (CONFIG_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${CONFIG_INSTALL_SUBDIR}")

# build shared libraries
set (BUILD_SHARED_LIBS ON)

# global property to hold the names of katydid library targets
set_property (GLOBAL PROPERTY ${PROJECT_NAME}_LIBRARIES)

# deal with the rpath settings (from http://www.cmake.org/Wiki/CMake_RPATH_handling)
# use (i.e. don't skip) the full RPATH for the build tree
#SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
#SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

#SET(CMAKE_INSTALL_RPATH "${LIB_INSTALL_DIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# add to the RPATH to be used when installing, but only if it's not a system directory
list (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${LIB_INSTALL_DIR}" isSystemDir)
if ("${isSystemDir}" STREQUAL "-1")
   list (APPEND CMAKE_INSTALL_RPATH "${LIB_INSTALL_DIR}")
endif ("${isSystemDir}" STREQUAL "-1")


##########
# MACROS #
##########

# This should be called immediately after setting the project name
macro (pbuilder_prepare_project VERSION_MAJOR VERSION_MINOR REVISION)
    # define the variables to describe the package (will go in the KatydidConfig.hh file)
    set (${PROJECT_NAME}_VERSION_MAJOR ${VERSION_MAJOR})
    set (${PROJECT_NAME}_VERSION_MINOR ${VERSION_MINOR})
    set (${PROJECT_NAME}_REVISION ${REVISION})
    set (${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_REVISION}")
    set (${PROJECT_NAME}_FULL_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_REVISION}")
    set (${PROJECT_NAME}_PACKAGE_NAME "${PROJECT_NAME}")
    set (${PROJECT_NAME}_PACKAGE_STRING "${PROJECT_NAME} ${${PROJECT_NAME}_FULL_VERSION}")
    
    # Configuration header file
    if (EXISTS ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.hh.in)
        configure_file (
            ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.hh.in
            ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.hh
        )
        # Add the binary tree to the search path for include files so that the config file is found during the build
        include_directories (${PROJECT_BINARY_DIR})
    endif (EXISTS ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.hh.in)
endmacro ()

macro (pbuilder_add_ext_libraries)
    list (APPEND EXTERNAL_LIBRARIES ${ARGN})
endmacro ()

macro (pbuilder_add_submodule_libraries)
    list (APPEND EXTERNAL_LIBRARIES ${ARGN})
    list (APPEND SUBMODULE_LIBRARIES ${ARGN})
endmacro ()

macro (pbuilder_install_libraries)
    install (TARGETS ${ARGN} EXPORT ${PROJECT_NAME}Targets DESTINATION ${LIB_INSTALL_DIR})
    #list (APPEND ${PROJECT_NAME}_LIBRARIES ${ARGN})
    set_property (GLOBAL APPEND PROPERTY ${PROJECT_NAME}_LIBRARIES ${ARGN})
    set_target_properties (${ARGN} PROPERTIES INSTALL_NAME_DIR ${LIB_INSTALL_DIR})
endmacro ()

macro (pbuilder_install_executables)
    install(TARGETS ${ARGN} EXPORT ${PROJECT_NAME}Targets DESTINATION ${BIN_INSTALL_DIR})
endmacro ()

macro (pbuilder_install_headers)
    install(FILES ${ARGN} DESTINATION ${INCLUDE_INSTALL_DIR})
endmacro ()

macro (pbuilder_install_header_dirs FILE_PATTERN)
    install(DIRECTORY ${ARGN} DESTINATION ${INCLUDE_INSTALL_DIR} FILES_MATCHING PATTERN "${FILE_PATTERN}")
endmacro ()

macro (pbuilder_install_config)
    install(FILES ${ARGN} DESTINATION ${CONFIG_INSTALL_DIR})
endmacro ()

macro (pbuilder_install_files DEST_DIR)
    install(FILES ${ARGN} DESTINATION ${DEST_DIR})
endmacro ()

# This should be called AFTER all subdirectories with libraries have been called, and all include directories added.
macro (pbuilder_install_config_files)
    # Configuration header file
    if (EXISTS ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.hh)
        # Install location for the configuration header
        pbuilder_install_headers (${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.hh)
    endif (EXISTS ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.hh)

    # CMake configuration file
    get_property(${PROJECT_NAME}_LIBRARIES GLOBAL PROPERTY ${PROJECT_NAME}_LIBRARIES)
    if (EXISTS ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in)
        # the awkwardness of the following four lines is because for some reason cmake wouldn't just go from .cmake.tmp to .cmake when I tested it.
        configure_file(${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}Config.cmake.tmp @ONLY)
        configure_file(${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}Config.cmake.tmp ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}Config.cmake.ppp @ONLY)
        file (REMOVE ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}Config.cmake.tmp)
        file (RENAME ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}Config.cmake.ppp ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}Config.cmake)
    endif (EXISTS ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in)
endmacro ()

macro (pbuilder_variables_for_parent)
    if (NOT ${PBUILDER_STANDALONE})
        get_property (LIBRARIES GLOBAL PROPERTY ${PROJECT_NAME}_LIBRARIES)
        set (${PROJECT_NAME}_LIBRARIES ${LIBRARIES} ${SUBMODULE_LIBRARIES} PARENT_SCOPE)
        set (${PROJECT_NAME}_LIBRARY_DIR ${LIB_INSTALL_DIR} PARENT_SCOPE)
        set (${PROJECT_NAME}_INCLUDE_DIR ${INCLUDE_INSTALL_DIR} PARENT_SCOPE)
    endif (NOT ${PBUILDER_STANDALONE})
endmacro ()

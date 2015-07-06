cmake_minimum_required( VERSION 2.8 )

#########
# setup #
#########

project( locust_mc )

list (APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

include( ${PROJECT_SOURCE_DIR}/PackageBuilder.cmake )
pbuilder_prepare_project( 1 0 0 )
include_directories( ${PROJECT_SOURCE_DIR} )

include( /usr/share/cmake-2.8/Modules/ExternalProject.cmake )
#include( ${CMAKE_MODULE_PATH}/ExternalProject.cmake )


#######################
# locust dependencies #
#######################


pbuilder_add_ext_libraries( "m" )

include_directories( ${PROJECT_SOURCE_DIR}/RapidJSON )

add_subdirectory( RandomLib )
include( ${CMAKE_BINARY_DIR}/RandomLib/randomlib-config.cmake )
include_directories( ${RandomLib_INCLUDE_DIRS} )
link_directories( ${RandomLib_LIBRARY_DIRS} )
pbuilder_add_ext_libraries( ${RandomLib_LIBRARIES} )
message( STATUS "randomlib library dirs: ${RandomLib_LIBRARY_DIRS}" )
# add to the RPATH to be used when installing, but only if it's not a system directory
list( FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${RandomLib_LIBRARY_DIRS}" isSystemDir )
message( STATUS "is system dir: ${isSystemDir}" )
if( "${isSystemDir}" STREQUAL "-1" )
   message( STATUS "appending to rpath" )
   list( APPEND CMAKE_INSTALL_RPATH ${RandomLib_LIBRARY_DIRS} )
endif( "${isSystemDir}" STREQUAL "-1" )

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
include_directories( ${PROJECT_SOURCE_DIR}/monarch/Source )
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



ExternalProject_Add( kassiopeia SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/kassiopeia )

pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaBindings.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaFields.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaGenerators.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaGeometry.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaInteractions.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaMath.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaModifiers.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaNavigators.so)
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaObjects.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaOperators.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaReaders.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaReaders.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaSimulation.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaTerminators.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaTrajectories.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaUtility.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaVisualization.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaWriters.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaUtility.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMCore.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMElectrostaticBoundaryIntegrals.so)
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMFileManipulation.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMGeoBagPlugin.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMHashGenerator.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMIntegratingFieldSolver.so)
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMIO.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMMath.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMStructuredASCII.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMSurfaces.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMVisualization.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMZHGenerator.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKEMZHSolver.so )

pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKommonCore.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaOperators.so )
pbuilder_add_ext_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/cbuild/kassiopeia-prefix/src/install/lib/libKassiopeiaReaders.so )

###########
# library #
###########

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

set (KASSIOPEIA_SOURCE_FILES
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSEvent.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSEventMessage.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSMainMessage.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRoot.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootElectricField.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootGenerator.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootMagneticField.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootSpace.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootSpaceInteraction.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootSpaceNavigator.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootStepModifier.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootSurfaceInteraction.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootSurfaceNavigator.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootTerminator.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootTrajectory.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRootWriter.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRun.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSRunMessage.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSSimulation.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSStep.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSStepMessage.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSTrack.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Source/KSTrackMessage.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kommon/Core/Utility/KMessage.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Source/KSCommand.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Source/KSCommandGroup.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Source/KSComponent.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Source/KSComponentGroup.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Source/KSDictionary.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Source/KSObject.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Source/KSObjectsMessage.cxx
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Source/KSToolbox.cxx
)

set (KASSIOPEIA_HEADER_FILES
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSEvent.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSEventMessage.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSMainMessage.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRoot.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootElectricField.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootGenerator.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootMagneticField.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootSpace.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootSpaceInteraction.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootSpaceNavigator.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootStepModifier.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootSurfaceInteraction.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootSurfaceNavigator.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootTerminator.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootTrajectory.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRootWriter.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRun.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSRunMessage.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSSimulation.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSStep.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSStepMessage.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSTrack.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Simulation/Include/KSTrackMessage.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kommon/Core/Utility/KMessage.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSCommandGroup.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSCommand.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSCommandMember.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSComponentGroup.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSComponent.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSComponentIntegral.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSCommandGroup.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSComponentMaximum.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSComponentMember.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSComponentMinimum.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSComponentTemplate.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSCore.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSDictionary.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSObject.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSObjectsMessage.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSOutput.h
${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/Objects/Include/KSToolbox.h
)


include_directories( ${KASSIOPEIA_INCLUDE_DIRECTORIES} )
#link_directories( ${KASSIOPEIA_INCLUDE_DIRECTORIES} )
#file(GLOB_RECURSE KASSIOPEIA_KEMFIELD ${PROJECT_SOURCE_DIR}/kassiopeia/Kassiopeia/KEMField/  *.cxx )
#add_library(Kassiopeia_KEMField ${KASSIOPEIA_KEMFIELD})
#add_library(Kassiopeia ${KASSIOPEIA_SOURCE_FILES})
#pbuilder_add_ext_libraries ( Kassiopeia )
#pbuilder_add_ext_libraries ( Kassiopeia_KEMField )
#pbuilder_install_headers( ${KASSIOPEIA_HEADER_FILES} )
#pbuilder_install_libraries( Kassiopeia )
#pbuilder_install_libraries( Kassiopeia_KEMField )


set( LOCUST_MC_HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCConfigurator.hh
#    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCConstants.hh
#    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCDefaultConfig.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCDestroyer.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCEggWriter.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCException.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCFactory.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCGenerator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCGeneratorToolbox.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCLogger.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCParam.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCParser.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCRunLengthCalculator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCSignal.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCSimulationController.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCSingleton.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCVisitor.hh
    
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCDigitizer.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCGaussianNoiseGenerator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCTestSignalGenerator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCBasebandTrackGenerator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCReceiverTransferFunctionGenerator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCTrappedElectronGenerator.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCKassiopeiaGenerator.hh
#    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCSingleAmpSim.hh
)

set( LOCUST_MC_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCConfigurator.cc
#    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCDefaultConfig.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCEggWriter.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCException.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCGenerator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCGeneratorToolbox.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCLogger.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCParam.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCParser.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCRunLengthCalculator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCSignal.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCSimulationController.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Core/LMCVisitor.cc

    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCDigitizer.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCGaussianNoiseGenerator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCTestSignalGenerator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCBasebandTrackGenerator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCReceiverTransferFunctionGenerator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCTrappedElectronGenerator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCKassiopeiaGenerator.cc
#    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Generators/LMCSingleAmpSim.cc
)

add_library( LocustMC ${LOCUST_MC_SOURCE_FILES} )
target_link_libraries( LocustMC ${EXTERNAL_LIBRARIES} )


pbuilder_install_headers( ${LOCUST_MC_HEADER_FILES} )
pbuilder_install_libraries( LocustMC )


###############
# executables #
###############


add_executable( LocustSim ${CMAKE_CURRENT_SOURCE_DIR}/Source/Applications/LocustSim.cc )
target_link_libraries( LocustSim LocustMC ${EXTERNAL_LIBRARIES} )
pbuilder_install_executables( LocustSim )


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


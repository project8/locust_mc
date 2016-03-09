if( locust_mc_ENABLE_EXECUTABLES )

    set( lib_dependencies 
        LocustMC
        )

IF (locust_mc_BUILD_WITH_KASSIOPEIA)
    set( lib_dependencies 
        LocustMC
KommonCore
KommonRoot
KommonGsl
KommonVtk
KGeoBagMath
KGeoBagMathLinearAlgebra
KGeoBagMathSpaceTree
KGeoBagCore
KGeoBagShapes
KGeoBagMesh
KGeoBagAxialMesh
KGeoBagDiscreteRotationalMesh
KGeoBagAppearance
KGeoBagMetrics
KGeoBagRandom
KGeoBagIntersection
KGeoBagDeformation
KGeoBagBasicVisualization
KGeoBagROOTVisualization
KGeoBagVtkVisualization
KGeoBagBindings
KEMCore
KEMHashGenerator
KEMStructuredASCII
KEMIO
KEMFileManipulation
KEMMath
KFMMathUtilities
KVectorMath
KEMSurfaces
KEMElectrostaticBoundaryIntegrals
KEMElectromagnets
KEMIntegratingFieldSolver
KEMZHGenerator
KEMZHSolver
KFMCore
KFMMath
KFMKernel
KFMTree
KFMElectrostatics
KFMInterfaceExtraction
KFMInterfaceFieldSolvers
KFMUtility
KEMVisualization
KEMGeoBagPlugin
KassiopeiaUtility
KassiopeiaMath
KassiopeiaObjects
KassiopeiaOperators
KassiopeiaFields
KassiopeiaGeometry
KassiopeiaGenerators
KassiopeiaTrajectories
KassiopeiaNavigators
KassiopeiaInteractions
KassiopeiaTerminators
KassiopeiaModifiers
KassiopeiaWriters
KassiopeiaReaders
KassiopeiaVisualization
KassiopeiaSimulation
KassiopeiaBindings
    )
ENDIF ( locust_mc_BUILD_WITH_KASSIOPEIA )


    
    set( programs
        LocustSim

    )
    
    pbuilder_executables( programs lib_dependencies )

endif( locust_mc_ENABLE_EXECUTABLES )

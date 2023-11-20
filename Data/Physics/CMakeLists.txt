# CMakeLists.txt for Locust_mc/Data/Physics

set( LOCUST_MC_DATA_FILES      
    KrOscillatorStrength.txt
    KrShakeParameters214.txt
    H2OscillatorStrength.txt
)

pbuilder_install_data( ${LOCUST_MC_DATA_FILES} )


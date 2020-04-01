# CMakeLists.txt for Locust_mc/Data

set( LOCUST_MC_DATA_FILES
    H2OscillatorStrength.txt
    KrOscillatorStrength.txt
    KrShakeParameters214.txt
)

pbuilder_install_data( ${LOCUST_MC_DATA_FILES} )

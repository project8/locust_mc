# CMakeLists.txt for Locust_mc/Data/ModeMaps

set( LOCUST_MC_DATA_FILES      
    PozarExport_Ex_RectangularTE011.txt
    PozarExport_Etheta_CylindricalTE011.txt
)

pbuilder_install_data( ${LOCUST_MC_DATA_FILES} )


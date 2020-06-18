# CMakeLists.txt for locust_mc/Config
# Author: N.S. Oblath

configure_file( JustKass.xml.in JustKass.xml )
configure_file( JustKassFieldMap.xml.in JustKassFieldMap.xml )
configure_file( LocustKass_Waveguide_Template.xml.in LocustKass_Waveguide_Template.xml )
configure_file( LocustKass_FreeSpace_Template.xml.in LocustKass_FreeSpace_Template.xml )
configure_file( LocustWaveguideTemplate.json.in LocustWaveguideTemplate.json )
configure_file( LocustFreeSpaceTemplate.json.in LocustFreeSpaceTemplate.json )
configure_file( LocustPlaneWaveTemplate.json.in LocustPlaneWaveTemplate.json )
configure_file( LocustFakeTrack.json.in LocustFakeTrack.json )
configure_file( LocustMagDipoleAntennaTemplate.json.in LocustMagDipoleAntennaTemplate.json )
configure_file( LocustTestSignal.json.in LocustTestSignal.json )
configure_file( plotfieldmap.c.in plotfieldmap.c )
configure_file( katydid.json.in katydid.json )


set (EXAMPLE_CONFIGFILES
    ${CMAKE_CURRENT_BINARY_DIR}/JustKass.xml
    ${CMAKE_CURRENT_BINARY_DIR}/JustKassFieldMap.xml
    ${CMAKE_CURRENT_BINARY_DIR}/LocustKass_Waveguide_Template.xml
    ${CMAKE_CURRENT_BINARY_DIR}/LocustKass_FreeSpace_Template.xml 
    WaveguideGeometry.xml
    FreeSpaceGeometry.xml
    ${CMAKE_CURRENT_BINARY_DIR}/LocustWaveguideTemplate.json    
    ${CMAKE_CURRENT_BINARY_DIR}/LocustFreeSpaceTemplate.json    
    ${CMAKE_CURRENT_BINARY_DIR}/LocustPlaneWaveTemplate.json    
    ${CMAKE_CURRENT_BINARY_DIR}/LocustFakeTrack.json    
    ${CMAKE_CURRENT_BINARY_DIR}/LocustMagDipoleAntennaTemplate.json    
    ${CMAKE_CURRENT_BINARY_DIR}/LocustTestSignal.json    
    ${CMAKE_CURRENT_BINARY_DIR}/katydid.json    
    ${CMAKE_CURRENT_BINARY_DIR}/plotfieldmap.c    
)

pbuilder_install_config (${EXAMPLE_CONFIGFILES})
